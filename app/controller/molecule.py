from typing import List, Optional
from sqlalchemy.orm import Session

import rdkit.Chem as Chem
from rdkit.Chem import Descriptors
from fastapi import UploadFile
import io
import re

from app.models.molecule import Molecule as MoleculeModel
from app.schemas.molecule import MoleculeCreate as SchemaMoleculeCreate
from app.schemas.molecule import MoleculeUpdate as SchemaMoleculeUpdate


def get_list_molecule(db: Session, skip: int = 0, limit: int = 100) -> List[MoleculeModel]:
    """get list of molecule from db"""
    return db.query(MoleculeModel).order_by(MoleculeModel.updated.desc()).offset(skip).limit(limit).all()


def generate_name(db: Session):
    count = db.query(MoleculeModel).filter(MoleculeModel.name.contains('molecule')).count()
    return f'molecule_{count+1}'


def create_molecule(db: Session, molecule: SchemaMoleculeCreate) -> MoleculeModel:
    """create molecule and save it to db"""
    properties = get_properties(molecule.smiles)
    db_molecule = MoleculeModel(
        smiles=molecule.smiles, name=molecule.name, description=molecule.description, 
        log_p=properties['log_p'], number_of_atoms=properties['number_of_atoms'], 
        molecular_weight=properties['molecular_weight']
    )
    try:
        db.add(db_molecule)
        db.commit()
    except Exception:
        db.rollback()
        
    db.refresh(db_molecule)
    return db_molecule


def get_properties(smiles):
    lig = Chem.MolFromSmiles(smiles)
    properties = {
        "number_of_atoms": lig.GetNumAtoms(),
        "molecular_weight": Descriptors.ExactMolWt(lig),
        "log_p": Descriptors.MolLogP(lig),
    }
    return properties


async def upload_molecule(db: Session, file: UploadFile):
    """read molecules in '.smi' file and save it to db"""
    molecules = []
    with file.file as f:
        for line in io.TextIOWrapper(f, encoding='utf-8'):
            space_regex = r"\s+"
            removed_space_line = re.sub(space_regex, ' ', line)
            line_split = removed_space_line.split(' ')
            line_split.pop()
            smile = line_split[-1]
            if len(line_split) > 1:
                name = removed_space_line.rsplit(' ', 2)[0]
            else:
                name = generate_name(db)
            molecules.append({
                'name': name,
                'smile': smile
            })
    
    for molecule in molecules:
        try:
            properties = get_properties(molecule['smile'])
            db_molecule = MoleculeModel(
                smiles=molecule['smile'], name=molecule['name'], log_p=properties['log_p'], 
                number_of_atoms=properties['number_of_atoms'], molecular_weight=properties['molecular_weight']
            )
            db.add(db_molecule)
            db.commit()
        except Exception:
            db.rollback()


def get_molecule(db: Session, molecule_id: int) -> Optional[MoleculeModel]:
    """get one molecule from db"""
    return db.query(MoleculeModel).filter(MoleculeModel.id == molecule_id).first()


def update_molecule(db: Session, molecule_id: int, molecule: SchemaMoleculeUpdate) -> Optional[MoleculeModel]:
    """update molecule and save it db"""
    db_molecule = db.query(MoleculeModel).filter(MoleculeModel.id == molecule_id).first()
    if db_molecule:
        for column, value in molecule.dict(exclude_unset=True).items():
            setattr(db_molecule, column, value)
        db.commit()
    return db_molecule


def delete_molecule(db: Session, molecule_id: int) -> Optional[MoleculeModel]:
    """delete molecule from db"""
    db_molecule = db.query(MoleculeModel).filter(MoleculeModel.id == molecule_id).first()
    if db_molecule:
        db.delete(db_molecule)
        db.commit()
    return db_molecule

