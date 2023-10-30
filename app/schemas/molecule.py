import datetime
from typing import Union

from pydantic import BaseModel


class MoleculeBase(BaseModel):
    """base schema for model molecule"""
    name: Union[str, None] = None
    description: Union[str, None] = None


class MoleculeCreate(MoleculeBase):
    """schema for create and update"""
    
    smiles: str

    pass


class MoleculeUpdate(MoleculeBase):
    """schema for create and update"""
    
    pass


class Molecule(MoleculeBase):
    """schema for list"""

    id: int
    smiles: str
    created: datetime.datetime
    updated: datetime.datetime

    class Config:
        orm_mode = True


class MoleculeDetail(MoleculeBase):
    """schema for detail"""

    id: int
    smiles: str
    log_p: float
    number_of_atoms: int = None
    molecular_weight: float = None
    created: datetime.datetime
    updated: datetime.datetime

    class Config:
        orm_mode = True
