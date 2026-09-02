import requests

ZOTERO_USER_ID = "3750207"      # numeric Zotero user ID, NOT your username
API_KEY = "xxx"

headers = {
    "Zotero-API-Key": API_KEY,
    "Content-Type": "application/json",
    "Accept": "application/json",
}

payload = [{"name": "DIM"}]

r = requests.post(
    f"https://api.zotero.org/users/{ZOTERO_USER_ID}/collections",
    headers=headers,
    json=payload,
    timeout=30,
)

print("status:", r.status_code)
print(r.text)
r.raise_for_status()
print(r.json())
