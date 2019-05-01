from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String
from sqlalchemy.types import Float
from sqlalchemy import ForeignKey
from sqlalchemy.orm import relationship
import numpy as np

#from zope.sqlalchemy import ZopeTransactionExtension
#DBSession = scoped_session(
#            sessionmaker(extension=ZopeTransactionExtension()))

Base = declarative_base()

class HuelsenbeckSimulationSet(Base):
    __tablename__ = 'huelsenbeck_sim_sets'
    id = Column(Integer, primary_key=True)
    rows = Column(Integer)
    cols = Column(Integer)
    theta = Column(Float)
    indelible_model = Column(String)
    genes = Column(Integer)
    m = Column(Integer)
    n = Column(Integer)
    a_max = Column(Float)
    b_max = Column(Float)
    nr_sims = Column(Integer)
    #huelsenbeck_simulations = relationship("HuelsenbeckSimulation", back_populates="simulation_set", lazy='joined')
    huelsenbeck_simulations = relationship("HuelsenbeckSimulation", back_populates="simulation_set")

class HuelsenbeckReconstructionSet(Base):
    __tablename__ = 'huelsenbeck_reconstruction_sets'
    id = Column(Integer, primary_key=True)
    sim_set_id = Column(Integer, ForeignKey('huelsenbeck_sim_sets.id'))
    method = Column(String)
    alignment_method = Column(String)
    kvd_method = Column(String)
    k = Column(String)
    simulation_set = relationship("HuelsenbeckSimulationSet", lazy='joined')
    #simulation_set = relationship("HuelsenbeckSimulationSet")
    huelsenbeck_reconstructions = relationship("HuelsenbeckReconstruction", back_populates="reconstruction_set", lazy='joined')

class HuelsenbeckSimulation(Base):
    __tablename__ = 'huelsenbeck_sims'
    id = Column(Integer, primary_key=True)
    sim_set_id = Column(Integer, ForeignKey('huelsenbeck_sim_sets.id'))
    sim_id = Column(Integer, ForeignKey('sims.id'))
    row = Column(Integer)
    col = Column(Integer)
    simulation = relationship("Simulation", lazy='joined')
    simulation_set = relationship("HuelsenbeckSimulationSet", lazy='joined', back_populates="huelsenbeck_simulations")

class HuelsenbeckReconstruction(Base):
    __tablename__ = 'huelsenbeck_reconstructions'
    id = Column(Integer, primary_key=True)
    reconstruction_set_id = Column(Integer, ForeignKey('huelsenbeck_reconstruction_sets.id'))
    reconstruction_id = Column(Integer, ForeignKey('reconstructions.id'))
    huelsenbeck_sim_id = Column(Integer, ForeignKey('huelsenbeck_sims.id'))
    huelsenbeck_simulation = relationship("HuelsenbeckSimulation", lazy='joined')
    reconstruction = relationship("Reconstruction", lazy='joined')
    reconstruction_set = relationship("HuelsenbeckReconstructionSet", back_populates="huelsenbeck_reconstructions", lazy='joined')

class Simulation(Base):
    __tablename__ = 'sims'
    id = Column(Integer, primary_key=True)
    tree = Column(String)
    theta = Column(Float)
    indelible_model = Column(String)
    genes = Column(Integer)
    m = Column(Integer)
    kmer_distance_matrices = relationship("SpeciesKmerDistanceMatrix", back_populates="simulation", lazy='joined')
    def __repr__(self):
        return "<Simulation(m='%d', genes='%d', t_a='%f', t_b='%f')>" % ( self.m, self.genes, self.t_a, self.t_b)

class SpeciesKmerDistanceMatrix(Base):
    __tablename__ = 'kmer_distance_matrices'
    id = Column(Integer, primary_key=True)
    sim_id = Column(Integer, ForeignKey('sims.id'))
    alignment_method = Column(String)
    method = Column(String)
    k = Column(Integer)
    D12 = Column(Float)
    D13 = Column(Float)
    D14 = Column(Float)
    D15 = Column(Float)
    D23 = Column(Float)
    D24 = Column(Float)
    D25 = Column(Float)
    D34 = Column(Float)
    D35 = Column(Float)
    D45 = Column(Float)
    simulation = relationship("Simulation", back_populates="kmer_distance_matrices", lazy='joined')
    def to_np(self):
        x = np.zeros((5,5))
        x[0,1]=x[1,0]=self.D12
        x[0,2]=x[2,0]=self.D13
        x[0,3]=x[3,0]=self.D14
        x[0,4]=x[4,0]=self.D15
        x[1,2]=x[2,1]=self.D23
        x[1,3]=x[3,1]=self.D24
        x[1,4]=x[4,1]=self.D25
        x[2,3]=x[3,2]=self.D34
        x[2,4]=x[4,2]=self.D35
        x[3,4]=x[4,3]=self.D45
        return(x)
    def to_dm(self):
        DistanceMatrix(['A','B','C','D','E']

def kmer_distance_matrix_from_dm(x, sim_id, method, k):
    return(SpeciesKmerDistanceMatrix(sim_id=sim_id,method=method,k=k,D12=x['A','B'],D13=x['A','C'],D14=x['A','D'],D15=x['A','E'],D23=x['B','C'],D24=x['B','D'],D25=x['B','E'],D34=x['C','D'],D35=x['C','E'],D45=x['D','E']))
def kmer_distance_matrix_from_np(x, sim_id, method, k):
    return(SpeciesKmerDistanceMatrix(sim_id=sim_id,method=method,k=k,D12=x[0,1],D13=x[0,2],D14=x[0,3],D15=x[0,4],D23=x[1,2],D24=x[1,3],D25=x[1,4],D34=x[2,3],D35=x[2,4],D45=x[3,4]))

class Reconstruction(Base):
    __tablename__ = 'reconstructions'
    id = Column(Integer, primary_key=True)
    sim_id = Column(Integer, ForeignKey('sims.id'))
    method = Column(String)
    alignment_method = Column(String)
    kvd_method = Column(String)
    k = Column(String)
    splits = Column(String)
    success = Column(Integer)
    simulation = relationship("Simulation", lazy='joined')

indelible_models = {
'JC':
"""
  [submodel]    JC
""",
'LAV0.01a':
"""
  [submodel]    JC
  [indelmodel]   LAV  1.1 100  //  Lavelette (a=1.1, M=100)
  [insertrate]   0.01       // relative to substitution rate of 1
  [deleterate]   0.01
""",
'LAV0.05a':
"""
  [submodel]    JC
  [indelmodel]   LAV  1.1 100 
  [insertrate]   0.05  
  [deleterate]   0.05
""",
'LAV0.1a':
"""
  [submodel]    JC
  [indelmodel]   LAV  1.1 100  
  [insertrate]   0.1       
  [deleterate]   0.1
""",
'LAV0.01b':
"""
  [submodel]    JC
  [indelmodel]   LAV  1.5 100  //  Lavelette (a=1.5, M=100)
  [insertrate]   0.01       // relative to substitution rate of 1
  [deleterate]   0.01
""",
'LAV0.05b':
"""
  [submodel]    JC
  [indelmodel]   LAV  1.5 100 
  [insertrate]   0.05  
  [deleterate]   0.05
""",
'LAV0.1b':
"""
  [submodel]    JC
  [indelmodel]   LAV  1.5 100  
  [insertrate]   0.1       
  [deleterate]   0.1
""",
'LAV0.01c':
"""
  [submodel]    JC
  [indelmodel]   LAV  1.8 100  //  Lavelette (a=1.8, M=100)
  [insertrate]   0.01       // relative to substitution rate of 1
  [deleterate]   0.01
""",
'LAV0.05c':
"""
  [submodel]    JC
  [indelmodel]   LAV  1.8 100 
  [insertrate]   0.05  
  [deleterate]   0.05
""",
'LAV0.1c':
"""
  [submodel]    JC
  [indelmodel]   LAV  1.8 100  
  [insertrate]   0.1       
  [deleterate]   0.1
""",
}
