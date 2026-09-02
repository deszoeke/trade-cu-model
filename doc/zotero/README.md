# Python JSON queries to the zotero web API.

Simon de Szoeke with assistance of GitHub copilot LLMs

Get a write permission Zotero API key by logging in to Zotero at 
[https://www.zotero.org/settings/security](https://www.zotero.org/settings/security).
Put the key in `./.env` (secret and gitignored):
```
ZOTERO_API_KEY = "my key here"
```

Create and activate an environment with requests and python-dotenv:
```
mamba create -n zotero requests python-dotenv
mamba activate zotero
```

Run `createcollection.py` _once_ to create a zotero collection DIM and get its key.

`query.py` queries first author, year pairs in my Zotero library and if found moves them to 
the collection.
