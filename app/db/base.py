from datetime import datetime

from sqlalchemy import DateTime, Integer
from sqlalchemy.ext.compiler import compiles
from sqlalchemy.orm import declarative_base
from sqlalchemy.sql import expression
from sqlalchemy.orm import Mapped, mapped_column


class utcnow(expression.FunctionElement):
    """declare type utcnow as DateTime"""

    type = DateTime()


@compiles(utcnow, "postgresql")
def pg_utcnow(element, compiler, **kw) -> str:
    """specific command for postgresql for getting utc"""
    return "TIMEZONE('utc', CURRENT_TIMESTAMP)"


Base = declarative_base()


class BaseModel(Base):
    """
    This class is the abstract model that inherited from other model
    """

    __abstract__ = True

    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    created: Mapped[datetime] = mapped_column(
        DateTime,
        default=datetime.utcnow,
        server_default=utcnow(),
        nullable=False,
    )
    updated: Mapped[datetime] = mapped_column(
        DateTime,
        default=datetime.utcnow,
        onupdate=datetime.utcnow,
        nullable=False,
        server_default=utcnow(),
        server_onupdate=utcnow(),
    )
