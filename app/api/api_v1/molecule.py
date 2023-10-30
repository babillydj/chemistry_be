from typing import List, Dict
from fastapi import APIRouter, Depends, HTTPException, UploadFile
from fastapi.responses import JSONResponse
from sqlalchemy.orm import Session

from app.controller import molecule as controller
from app.db.session import get_db
from app.schemas import molecule as schemas
from app.models.molecule import Molecule as MoleculeModel

router =  APIRouter(
    prefix="/api/v1",
)


@router.get("/molecule", tags=["molecule"], response_model=List[schemas.Molecule])
def molecule_list(skip: int = 0, limit: int = None, sort_by: str = 'updated', db: Session = Depends(get_db)) -> List[MoleculeModel]:
    molecules = controller.get_list_molecule(db, skip=skip, limit=limit, sort_by=sort_by)
    return molecules


@router.post("/molecule", tags=["molecule"])
def molecule_create(molecule: schemas.MoleculeCreate, db: Session = Depends(get_db)) -> JSONResponse:
    controller.create_molecule(db=db, molecule=molecule)
    return JSONResponse(content={"message": "molecule created"})


@router.get("/molecule/{molecule_id}", tags=["molecule"], response_model=schemas.MoleculeDetail)
def molecule_detail(molecule_id: int, db: Session = Depends(get_db)) -> MoleculeModel:
    db_molecule = controller.get_molecule(db, molecule_id=molecule_id)
    if db_molecule is None:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return db_molecule


@router.put("/molecule/{molecule_id}", tags=["molecule"])
def molecule_update(
    molecule_id: int, molecule: schemas.MoleculeUpdate, db: Session = Depends(get_db)
) -> JSONResponse:
    db_molecule = controller.update_molecule(db, molecule=molecule, molecule_id=molecule_id)
    if db_molecule is None:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return JSONResponse(content={"message": "molecule updated"})


@router.delete("/molecule/{molecule_id}", tags=["molecule"])
def molecule_delete(
    molecule_id: int, db: Session = Depends(get_db)
) -> JSONResponse:
    db_molecule = controller.delete_molecule(db, molecule_id=molecule_id)
    if db_molecule is None:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return JSONResponse(content={"message": "molecule deleted"})


@router.post("/upload-smile/")
async def smile_upload(file: UploadFile, db: Session = Depends(get_db)):
    ext_type = file.filename.split('.')
    if ext_type[-1] != 'smi':
        raise HTTPException(status_code=400, detail="File not supported")
    await controller.upload_molecule(db, file=file)
    return JSONResponse(content={"message": "molecule uploaded"})
