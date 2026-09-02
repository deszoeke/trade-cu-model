# Python JSON queries to the zotero web API.

Simon de Szoeke with assistance of GitHub copilot LLMs

Get a write permission Zotero API key and put it in .env (secret and gitignored).

Run `createcollection.py` _once_ to create a zotero collection DIM and get its key.

`query.py` queries first author, year pairs in my Zotero library and if found moves them to 
the collection.
