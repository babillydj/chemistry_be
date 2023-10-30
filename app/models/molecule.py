from sqlalchemy import Integer, Float, String
from sqlalchemy.orm import Mapped, mapped_column

from app.db.base import BaseModel


class Molecule(BaseModel):
    """represent data model molecule"""

    __tablename__ = "molecules"
    smiles: Mapped[str] = mapped_column(String, unique=True)
    name: Mapped[str] =  mapped_column(String)
    description: Mapped[str] =  mapped_column(String, nullable=True)
    log_p: Mapped[float] = mapped_column(Float, default=0)
    number_of_atoms: Mapped[int] = mapped_column(Integer, default=0)
    molecular_weight: Mapped[float] = mapped_column(Float, default=0)
