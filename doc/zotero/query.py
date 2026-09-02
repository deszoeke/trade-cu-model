import os
import requests
from dotenv import load_dotenv

USER_ID = "3750207"
load_dotenv()
API_KEY = os.environ["ZOTERO_API_KEY"]
headers = {
    "Zotero-API-Key": API_KEY,
    "Content-Type": "application/json",
    "Accept": "application/json",
}

# reuse the already-created collection; do not create a new one
collection_key = "GHQFTWV4"
print("collection key:", collection_key)

pairs = [
    ("Albright", 2023),
    ("Bony", 2005),
    ("Bony", 2006),
    ("Bony", 2021),
    ("Bony", 2022),
    ("Bretherton", 2013),
    ("Brient", 2016),
    ("Coakley", 1975),
    ("George", 2022),
    ("Held", 2006),
    ("Jacob", 2020),
    ("Jansson", 2023),
    ("Jin", 2026),
    ("Lorenz", 2010),
    ("Mellor", 1977),
    ("Minnis", 2008),
    ("Myers", 2014),
    ("Nuijens", 2013),
    ("Nuijens", 2015),
    ("Richter", 2008),
    ("Sherwood", 2014),
    ("Siebesma", 2003),
    ("Sommeria", 1977),
    ("Stephan", 2021),
    ("Stevens", 2019),
    ("Stevens", 2021),
    ("Sun", 2022),
    ("Tornow", 2020),
    ("Tselioudis", 2021),
    ("Vecchi", 2007),
    ("Vial", 2017),
    ("Vogel", 2020),
    ("Vogel", 2022),
    ("Webb", 2012),
    ("Webb", 2013),
    ("Webb", 2017),
    ("Zelinka", 2012),
    ("Zelinka", 2023),
    ("Zelinka", 2025),
]

seen = set()
added = 0
not_found = []

# the API's q param is a plain phrase search (no creator:/year: field syntax,
# that only works in the desktop app's Advanced Search UI), so search by
# last name and filter the results by year and creator ourselves.
for author, year in pairs:
    r = requests.get(
        f"https://api.zotero.org/users/{USER_ID}/items",
        headers=headers,
        params={"format": "json", "limit": 50, "q": author, "qmode": "titleCreatorYear"},
        timeout=30,
    )
    r.raise_for_status()
    hits = r.json()

    matched = []
    for item in hits:
        data = item.get("data", {})
        if str(year) not in (data.get("date") or ""):
            continue
        creators = data.get("creators", [])
        if not creators:
            continue
        first = creators[0]
        first_name = first.get("lastName") or first.get("name") or ""
        if author.lower() in first_name.lower():
            matched.append(item)

    if not matched:
        not_found.append((author, year))
        continue

    for item in matched:
        item_key = item.get("key")
        if not item_key or item_key in seen:
            continue
        seen.add(item_key)

        data = item.get("data", {})
        collections = list(data.get("collections", []))
        if collection_key in collections:
            print(f"already in collection {author} {year}: {item_key}")
            added += 1
            continue
        collections.append(collection_key)

        # membership is set via the item's own collections field, there is
        # no /collections/{key}/items endpoint to POST to
        add = requests.patch(
            f"https://api.zotero.org/users/{USER_ID}/items/{item_key}",
            headers={**headers, "If-Unmodified-Since-Version": str(item.get("version"))},
            json={"collections": collections},
            timeout=30,
        )
        if add.status_code >= 400:
            print(f"FAILED {author} {year}: {item_key} -> {add.status_code} {add.text}")
            continue
        print(f"added {author} {year}: {item_key}")
        added += 1

print("total added:", added)
if not_found:
    print("no match found for:", not_found)


