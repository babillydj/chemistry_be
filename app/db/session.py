from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

engine = create_engine("postgresql+psycopg2://fa_user:fa_password@db/chemistry_be")

SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)


# Dependency
def get_db():
    """dependency db session"""
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()
