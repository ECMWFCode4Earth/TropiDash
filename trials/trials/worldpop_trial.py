# (failed) trials on downloading tiff files by API form worldpop database
# trials using json
import requests
import json

pop = requests.get("https://www.worldpop.org/rest/data/pop")
jpop = json.loads(pop.content)
jpop["data"]

# pop = requests.get("https://www.worldpop.org/rest/data/pop/wpgp")
# jpop = json.loads(pop.content)
# jpop["data"][0]
requests.get("https://www.worldpop.org/rest/data/pop_density/pd_ic_1km").json()

requests.get("https://www.worldpop.org/rest/data/pop").json()

requests.get("https://www.worldpop.org/rest/data/pop/wpgp1km/data").json()


#to try
requests.get("https://www.worldpop.org/rest/data/pop/wpgp1km?id=24777").json()
