import os 
from dotenv import load_dotenv

from sqlalchemy import (
    create_engine, Column, Integer, String, ForeignKey, Text
)
from sqlalchemy.orm import (
    declarative_base, sessionmaker, relationship
)

load_dotenv()
SUPABASE_DB_URL = os.getenv("SUPABASE_DB_URL")
print(f"SUPABASE_DB_URL: {SUPABASE_DB_URL}")
if not SUPABASE_DB_URL:
    raise RuntimeError("SUPABASE_DB_URL not set in environment")

engine = create_engine(
    SUPABASE_DB_URL,
    echo=False,       
    future=True      
)
SessionLocal = sessionmaker(
    bind=engine,
    autoflush=False,
    autocommit=False,
    expire_on_commit=False
)

Base = declarative_base()

class Sample(Base):
    __tablename__ = "samples"
    id       = Column(Integer, primary_key=True, index=True)
    sample_id = Column(String, unique=True, index=True)
    description = Column(Text, nullable=True)
    sequence    = Column(Text, nullable=False)
    mutations = relationship("Mutation", back_populates="sample")

class Mutation(Base):
    __tablename__ = "mutations"
    id            = Column(Integer, primary_key=True)
    sample_id     = Column(Integer, ForeignKey("samples.id"))
    reference_id  = Column(Integer, ForeignKey("reference_sequences.id"))  
    position      = Column(Integer)
    mutation      = Column(String)
    drug_class    = Column(String)
    interpretation= Column(Text)
    sample        = relationship("Sample", back_populates="mutations")
    reference     = relationship("ReferenceSequence")

class MutationRule(Base):
    __tablename__ = "mutation_rules"
    id          = Column(Integer, primary_key=True)
    position    = Column(Integer, nullable=False)
    ref_aa      = Column(String, nullable=False)
    alt_aa      = Column(String, nullable=False)
    drug_class  = Column(String)
    description = Column(Text)

class ReferenceSequence(Base):
    __tablename__ = "reference_sequences"
    id          = Column(Integer, primary_key=True)
    name        = Column(String, unique=True)
    description = Column(Text)
    sequence    = Column(Text)
    mutations   = relationship("Mutation", back_populates="reference")

def init_db():
    Base.metadata.create_all(bind=engine)

def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()
