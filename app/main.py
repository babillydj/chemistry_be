from typing import Dict
from fastapi import FastAPI

from app.api.api_v1 import molecule
from app.db.base import Base
from app.db.session import engine

from starlette.middleware.cors import CORSMiddleware

Base.metadata.create_all(bind=engine)

app = FastAPI()
app.include_router(molecule.router)

origins=['*']
app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/")
def read_root() -> Dict[str, str]:
    """open root"""
    return {"Hello": "World"}
