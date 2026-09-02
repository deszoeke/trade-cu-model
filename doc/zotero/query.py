import os
import requests

USER_ID = "3750207"
API_KEY = os.environ["ZOTERO_API_KEY"]  # never hardcode the key here

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

for author, year in pairs:
    q = f'creator:"{author}" year:{year}'
    r = requests.get(
        f"https://api.zotero.org/users/{USER_ID}/items",
        headers=headers,
        params={"format": "json", "limit": 20, "q": q, "qmode": "everything"},
        timeout=30,
    )
    r.raise_for_status()
    hits = r.json()

    if not hits:
        not_found.append((author, year))
        continue

    for item in hits:
        item_key = item.get("key")
        if not item_key or item_key in seen:
            continue
        seen.add(item_key)

        add = requests.post(
            f"https://api.zotero.org/users/{USER_ID}/collections/{collection_key}/items",
            headers=headers,
            json=[{"key": item_key}],
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


