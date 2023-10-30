from fastapi.testclient import TestClient
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from app.db.base import Base
from app.db.session import get_db
from app.main import app

from app.models.molecule import Molecule as MoleculeModel

SQLALCHEMY_DATABASE_URL = "sqlite:///./test.db"

engine = create_engine(
    SQLALCHEMY_DATABASE_URL, connect_args={"check_same_thread": False}
)
TestingSessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)


Base.metadata.create_all(bind=engine)


def override_get_db():
    db = None
    try:
        db = TestingSessionLocal()
        yield db
    finally:
        db.close()


app.dependency_overrides[get_db] = override_get_db

client = TestClient(app)
    

def test_create_item():
    response = client.post(
        "api/v1/molecule",
        json={"name": "Ethanol", "smiles": "CCO"},
    )
    assert response.status_code == 200
    assert response.json() == {
        "message": "molecule created"
    }


def test_get_list_molecule():
    response = client.get("api/v1/molecule")
    assert response.status_code == 200
    data = response.json()
    assert len(data) == 1
    assert data[0]["id"] == 1
    assert data[0]["name"] == "Ethanol"
    assert data[0]["smiles"] == "CCO"


def test_get_molecule():
    response = client.get("api/v1/molecule/1")
    assert response.status_code == 200
    data = response.json()
    assert data["id"] == 1
    assert data["name"] == "Ethanol"
    assert data["smiles"] == "CCO"


def test_update_molecule():
    response = client.put(
        "api/v1/molecule/1",
        json={"name": "Ethanol v1"},
    )
    assert response.status_code == 200
    assert response.json() == {
        "message": "molecule updated"
    }

    response = client.get("api/v1/molecule/1")
    assert response.status_code == 200
    data = response.json()
    assert data["id"] == 1
    assert data["name"] == "Ethanol v1"


def test_get_bad_molecule():
    response = client.get("api/v1/molecule/0")
    assert response.status_code == 404
    assert response.json() == {
        "detail": "Molecule not found"
    }


def test_delete_molecule():
    response = client.delete("api/v1/molecule/1")
    assert response.status_code == 200
    assert response.json() == {
        "message": "molecule deleted"
    }
    response = client.get("api/v1/molecule/1")
    assert response.status_code == 404
    

def test_upload_file():
    test_file = 'test.smi'
    files = {'file': ('test.smi', open(test_file, 'rb'))}
    response = client.post(
        "api/v1/upload-smile",
        files=files
    )
    assert response.status_code == 200
    assert response.json() == {
        "message": "molecule uploaded"
    }
    
    response = client.get("api/v1/molecule")
    assert response.status_code == 200
    data = response.json()
    assert len(data) == 5
    
def test_delete_all():
    db = TestingSessionLocal()
    db.query(MoleculeModel).delete()
    db.commit()
