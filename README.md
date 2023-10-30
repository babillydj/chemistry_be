# Description
this is a simple api FastAPI project,   


# Installation
Make sure you have the latest Docker and Docker Compose installed and active
- copy .env.dev to .env
- run `docker-compose build`
- run `docker-compose up`

It should be running at *localhost:8000/* now,     
you can open *localhost:8000/redoc/* or *localhost:8000/molecule/* for trial   


run test:   
- run `docker exec -it chemistry_be-server-1 sh` (container name)
- run `pytest test/molecule.py`
