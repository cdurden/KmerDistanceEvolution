from sqlalchemy import engine_from_config
from sqlalchemy.event import listen

from sqlalchemy.orm import scoped_session,sessionmaker,subqueryload,lazyload,Load
from sqlalchemy import and_
from sqlalchemy.ext.declarative import declarative_base
#Session = sessionmaker(bind=engine)
#session = Session()

from zope.sqlalchemy import ZopeTransactionExtension
DBSession = scoped_session(
                    sessionmaker(extension=ZopeTransactionExtension()))
from sqlalchemy import create_engine
import os

from KmerCoalescent.DB import *

engine = create_engine('sqlite:///db.sql', echo=True, convert_unicode=True)
from sqlalchemy.orm import sessionmaker
Session = sessionmaker(bind=engine)
session = Session()
Base.metadata.create_all(engine)

def huelsenbeck_reconstruction_set(reconstruction_set_id):
    import numpy as np
    reconstruction_set = session.query(HuelsenbeckReconstructionSet).get(reconstruction_set_id)
    return(reconstruction_set)

reconstruction_set = huelsenbeck_reconstruction_set(109)
for huel_rec in rec_set.huelsenbeck_reconstructions:
    print(huel_rec.reconstruction.success)


a_max = reconstruction_set.simulation_set.a_max
b_max = reconstruction_set.simulation_set.b_max
nr_rows = reconstruction_set.simulation_set.rows
nr_cols = reconstruction_set.simulation_set.cols

# Name the output file

# Read input data from FITS file
n = np.zeros((nr_rows, nr_cols))
d = np.zeros((nr_rows, nr_cols))
for huelsenbeck_reconstruction in reconstruction_set.huelsenbeck_reconstructions:
    row = huelsenbeck_reconstruction.huelsenbeck_simulation.row
    col = huelsenbeck_reconstruction.huelsenbeck_simulation.col
    d[row,col] += 1
    if bool(huelsenbeck_reconstruction.reconstruction.success):
        n[row,col] += 1
print(n)
print(d)
z = n/d

sim = reconstruction_set.simulation_set.huelsenbeck_simulations[0].simulation
for dm in sim.kmer_distance_matrices:
    if dm.k==1:
        print dm.distance_formula
        print dm.to_dm()

sim_set = reconstruction_set.simulation_set
huel_rec = reconstruction_set.huelsenbeck_reconstructions[0]
huel_sim = huel_rec.huelsenbeck_simulation
rec = huel_rec.reconstruction
sim = rec.simulation
for dm in sim.kmer_distance_matrices:
    if dm.k==1:
        print dm.distance_formula
        print dm.to_dm()



def huelsenbeck(request):
    try:
        methods = [request.matchdict['method']]
        distance_formulas = [request.matchdict['distance_formula']]
        k = [request.matchdict['k']]
    except KeyError:
        methods=request.GET.getall('method')
        distance_formulas=request.GET.getall('distance_formula')
        k = request.GET.getall('k')
    try:
        n = [request.matchdict['n']]
        m = [request.matchdict['m']]
        indelible_models = [request.matchdict['indelible_model']]
        genes = [request.matchdict['genes']]
        theta = [request.matchdict['theta']]
    except KeyError:
        n = request.GET.getall('n')
        m = request.GET.getall('m')
        indelible_models=request.GET.getall('indelible_model')
        theta = request.GET.getall('theta')
        genes = request.GET.getall('genes')
    print(distance_formulas)
    print(k)
    if len(k)>0:
        k = [",".join(sorted(k))]
    else:
        k = []
    filter_group = [col.in_(vals) for col,vals in [(HuelsenbeckReconstructionSet.method, methods),
                                                   (HuelsenbeckReconstructionSet.distance_formula, distance_formulas),
                                                   (HuelsenbeckReconstructionSet.k, k),
                                                   (HuelsenbeckSimulationSet.genes, genes),
                                                   (HuelsenbeckSimulationSet.indelible_model, indelible_models),
                                                   (HuelsenbeckSimulationSet.n, n),
                                                   (HuelsenbeckSimulationSet.theta, theta),
                                                   (HuelsenbeckSimulationSet.m, m)] if len(vals)>0]

    reconstruction_set_macro = get_renderer("templates/huelsenbeck_reconstruction_set.pt").implementation().macros
    huelsenbeck_reconstruction_sets = DBSession.query(HuelsenbeckReconstructionSet).\
            join(HuelsenbeckReconstructionSet.simulation_set, aliased=True).\
            options(subqueryload(HuelsenbeckReconstructionSet.simulation_set)).\
            filter(and_(*filter_group)).all()
    huelsenbeck_reconstruction_sets.sort(key=lambda huelsenbeck_reconstruction_set: sum(huelsenbeck_reconstruction.reconstruction.success for huelsenbeck_reconstruction in huelsenbeck_reconstruction_set.huelsenbeck_reconstructions),reverse=True)
    return({'methods': methods, 'distance_formulas': distance_formulas, 'k': ",".join(k), 'm': m, 'genes': genes, 'huelsenbeck_reconstruction_sets': huelsenbeck_reconstruction_sets, 'results': [], 'reconstruction_set_macro': reconstruction_set_macro, 'np': np})

    import os
    engine = engine_from_config(settings, 'sims.')
    listen(engine, "connect", do_connect)
    listen(engine, "begin", do_begin)
    DBSession.configure(bind=engine)
    Base.metadata.bind = engine
    Base.metadata.create_all(engine)
