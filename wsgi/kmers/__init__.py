from wsgiref.simple_server import make_server
from pyramid.config import Configurator
from pyramid.response import Response
from pyramid.httpexceptions import HTTPFound
from pyramid.renderers import render_to_response
from pyramid.renderers import get_renderer

#from sqlalchemy import create_engine
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
#Base = declarative_base()

from KmerDistanceEvolution.Data import *

def do_connect(dbapi_connection, connection_record):
    # disable pysqlite's emitting of the BEGIN statement entirely.
    # also stops it from emitting COMMIT before any DDL.
    dbapi_connection.isolation_level = None

def do_begin(conn):
    # emit our own BEGIN
    conn.execute("BEGIN")

def nonesorter(a):
    if not a:
        return ""
    return a

def huelsenbeck_svg(request):
    print("running huelsenbeck_svg view code")
    import numpy as np
    import tkinter
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import xticks,yticks
    import matplotlib.cm as cm
    from mpl_toolkits.axes_grid1 import ImageGrid
    from math import log
    import os
    print('svg dir')
    svg_dir = os.path.join(request.registry.settings['static_dir'],'tmp')
    print(svg_dir)
    import tempfile
    plt.ioff()
    cmap = cm.jet
    tree_estimate_set_id = request.matchdict['tree_estimate_set_id']
#    method=request.matchdict['method']
#    m=request.matchdict['m']
#    genes=request.matchdict['genes']
#    kvd_method=request.matchdict['kvd_method']
#    k=tuple(request.matchdict['k'].split("_"))

    tree_estimate_set = DBSession.query(HuelsenbeckTreeEstimateSet).get(tree_estimate_set_id)
    from sqlalchemy.orm import joinedload
    from sqlalchemy.orm import defaultload
    a_max = tree_estimate_set.simulation_set.a_max
    b_max = tree_estimate_set.simulation_set.b_max
    nr_rows = tree_estimate_set.simulation_set.rows
    nr_cols = tree_estimate_set.simulation_set.cols

    # set of arrays to store cell data
    n = np.zeros((nr_rows, nr_cols))
    d = np.zeros((nr_rows, nr_cols))
    # Method 1
#    for row in range(nr_rows):
#        print(row)
#        for col in range(nr_cols):
#            print(col)
#            huelsenbeck_tree_estimates = DBSession.query(HuelsenbeckTreeEstimate).filter(tree_estimate_set==tree_estimate_set).filter(row==row).filter(col==col).all()
#            n[row,col] = sum([huelsenbeck_tree_estimate.success for huelsenbeck_tree_estimate in huelsenbeck_tree_estimates])
#            d[row,col] = len(huelsenbeck_tree_estimates)
#

    # Name the output file


    # Read input data from FITS file
    for huelsenbeck_tree_estimate in tree_estimate_set.huelsenbeck_tree_estimates:
        row = huelsenbeck_tree_estimate.huelsenbeck_simulation.row
        col = huelsenbeck_tree_estimate.huelsenbeck_simulation.col
        d[row,col] += 1
        if bool(huelsenbeck_tree_estimate.tree_estimate.success):
            n[row,col] += 1
    print(n)
    print(d)
    z = n/d

    print(z)
    # Define x and y labels
    x = np.linspace(0,a_max,z.shape[1])
    y = np.linspace(0,b_max,z.shape[0])

    # Set colomap
    cmap = cm.jet

    # Setup contour levels
    levels = np.arange(0,1.1,0.1)

    # Setup figure and colorbar
    fig = plt.figure(figsize=(6,6))
    #ax = fig.add_subplot(111)
    ax = fig.add_axes((.1,.3,.8,.6))
    #fig = plt.figure(figsize=(6,8))
    #ax = plt.add_axes((.1,.3,.8,.6))

    #plt.locator_params(axis='y', tight=True, nbins=nr_rows-1)
    #plt.locator_params(axis='x', tight=True, nbins=nr_cols-1)


    #im = ax.contourf(x, y, z, levels=levels, antialiased=False, cmap=cmap)
    im = plt.imshow(z, origin='lower', clim = (0.0,1.0))
    ax.set_xlabel('a')
    ax.set_ylabel('b')
    avals = [(col+1)*(a_max)/nr_cols for col in range(nr_cols)]
    #np.linspace(0,a_max,nr_cols)
    bvals = [(row+1)*(b_max)/nr_rows for row in range(nr_rows)]
    #np.linspace(0,b_max,nr_rows)
    xticks(range(nr_cols), ["{:0.2f}".format(aval) for aval in list(avals)])
    yticks(range(nr_rows), ["{:0.2f}".format(bval) for bval in list(bvals)])
    #xticks(xticks()[0][1:-1], [str(t*a_max/float(nr_cols-1)) for t in xticks()[0][1:-1]])
    #yticks(yticks()[0][1:-1], reversed([str(t*b_max/float(nr_rows-1)) for t in yticks()[0][1:-1]]))

    im.set_cmap(cmap)

    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(im, cax=cax)

    #here = os.path.dirname(os.path.abspath(__file__))
    #fd, path = tempfile.mkstemp(suffix='.svg',dir=os.path.join(here,'static','tmp'))
    svg_dir = os.path.join(request.registry.settings['static_dir'],'tmp')
    print(svg_dir)
    fd, path = tempfile.mkstemp(suffix='.svg',dir=svg_dir)
    plt.savefig(path)
    plt.close()
    #return(svg_dir)
    return HTTPFound(request.static_path('kmers:'+os.path.join('static','tmp',os.path.split(path)[1])))

def results(request):
    sim_parameters = ['n','indelible_model','genes','m','theta']
    tree_estimate_parameters = ['k','distance_formula','alignment_method','method']
    simulation_sets = DBSession.query(HuelsenbeckSimulationSet).options(Load(HuelsenbeckSimulationSet).lazyload('*')).all()
    tree_estimate_sets = DBSession.query(HuelsenbeckTreeEstimateSet).options(Load(HuelsenbeckTreeEstimateSet).lazyload('*')).all()
    sim_parameters_dict = dict([(parameter,set([getattr(simulation_set,parameter) for simulation_set in simulation_sets])) for parameter in sim_parameters])
    tree_estimate_parameters_dict = dict([(parameter,set([getattr(tree_estimate_set,parameter) for tree_estimate_set in tree_estimate_sets])) for parameter in tree_estimate_parameters])
    import itertools
    sim_parameter_combos = list(itertools.product(*[sim_parameters_dict[sim_parameter] for sim_parameter in sim_parameters]))
    tree_estimate_parameter_combos = list(itertools.product(*[tree_estimate_parameters_dict[tree_estimate_parameter] for tree_estimate_parameter in tree_estimate_parameters]))
    sim_parameter_combos.sort(key=lambda x: x if isinstance(x, str) else str(x))
    tree_estimate_parameter_combos.sort(key=lambda x: x if isinstance(x, str) else str(x))
    results = np.empty((len(sim_parameter_combos),len(tree_estimate_parameter_combos)),dtype=tuple)
    nonempty_cols = list()
    for row,sim_parameter_combo in enumerate(sim_parameter_combos):
        for col,tree_estimate_parameter_combo in enumerate(tree_estimate_parameter_combos):
            q = DBSession.query(HuelsenbeckTreeEstimateSet).\
                    join(HuelsenbeckTreeEstimateSet.simulation_set, aliased=True).\
                    options(lazyload(HuelsenbeckTreeEstimateSet.huelsenbeck_tree_estimates))
            for param, value in zip(sim_parameters,sim_parameter_combo):
                    q = q.filter(getattr(HuelsenbeckSimulationSet, param)==value)
            for param, value in zip(tree_estimate_parameters,tree_estimate_parameter_combo):
                    q = q.filter(getattr(HuelsenbeckTreeEstimateSet, param)==value)
            results[row,col]= tuple([tree_estimate_set.id for tree_estimate_set in q.all()])
            if len(results[row,col])>0:
                nonempty_cols.append(col)
    nonempty_cols = list(set(nonempty_cols))
    results = results[:,nonempty_cols]
    tree_estimate_parameter_combos = [tree_estimate_parameter_combos[col] for col in nonempty_cols]
    return({'results': results, 'sim_parameters':sim_parameters, 'tree_estimate_parameters': tree_estimate_parameters, 'sim_parameter_combos':sim_parameter_combos, 'tree_estimate_parameter_combos':tree_estimate_parameter_combos})


def huelsenbeck_tree_estimate_set(request):
    import numpy as np

    tree_estimate_set_id = request.matchdict['tree_estimate_set_id']

    tree_estimate_set = DBSession.query(HuelsenbeckTreeEstimateSet).get(tree_estimate_set_id)
    return({'huelsenbeck_tree_estimate_set': tree_estimate_set, 'np': np})

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
    filter_group = [col.in_(vals) for col,vals in [(HuelsenbeckTreeEstimateSet.method, methods),
                                                   (HuelsenbeckTreeEstimateSet.distance_formula, distance_formulas),
                                                   (HuelsenbeckTreeEstimateSet.k, k),
                                                   (HuelsenbeckSimulationSet.genes, genes),
                                                   (HuelsenbeckSimulationSet.indelible_model, indelible_models),
                                                   (HuelsenbeckSimulationSet.n, n),
                                                   (HuelsenbeckSimulationSet.theta, theta),
                                                   (HuelsenbeckSimulationSet.m, m)] if len(vals)>0]

    tree_estimate_set_macro = get_renderer("templates/huelsenbeck_tree_estimate_set.pt").implementation().macros
    huelsenbeck_tree_estimate_sets = DBSession.query(HuelsenbeckTreeEstimateSet).\
            join(HuelsenbeckTreeEstimateSet.simulation_set, aliased=True).\
            options(subqueryload(HuelsenbeckTreeEstimateSet.simulation_set)).\
            filter(and_(*filter_group)).all()
    #huelsenbeck_tree_estimate_sets.sort(key=lambda huelsenbeck_tree_estimate_set: sum(huelsenbeck_tree_estimate.tree_estimate.success for huelsenbeck_tree_estimate in huelsenbeck_tree_estimate_set.huelsenbeck_tree_estimates),reverse=True)
    return({'methods': methods, 'distance_formulas': distance_formulas, 'k': ",".join(k), 'm': m, 'genes': genes, 'huelsenbeck_tree_estimate_sets': huelsenbeck_tree_estimate_sets, 'results': [], 'tree_estimate_set_macro': tree_estimate_set_macro, 'np': np})

if __name__ == '__main__':
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    config = Configurator()
    config.include('pyramid_chameleon')
    #config.add_route('huelsenbeck', '/huelsenbeck/{method}/{distance_formula}/{genes}/{m}/{k}')
    config.add_route('huelsenbeck_svg', '/huelsenbeck_svg/{tree_estimate_set_id}.svg')
    config.add_view(huelsenbeck_svg, route_name='huelsenbeck_svg', renderer=os.path.join(here,'templates','huelsenbeck_tree_estimate_set.pt'))
    config.add_route('huelsenbeck_tree_estimate_set', '/huelsenbeck_tree_estimate_set/{tree_estimate_set_id}')
    config.add_view(huelsenbeck_tree_estimate_set, route_name='huelsenbeck_tree_estimate_set', renderer=os.path.join(here,'templates','huelsenbeck_tree_estimate_set.pt'))
    config.add_route('huelsenbeck', '/huelsenbeck')
    config.add_view(huelsenbeck, route_name='huelsenbeck', renderer=os.path.join(here,'templates','huelsenbeck.pt'))
    config.add_route('results', '/results')
    config.add_view(results, route_name='results', renderer=os.path.join(here,'templates','results.pt'))

    app = config.make_wsgi_app()
    server = make_server('0.0.0.0', 8080, app)
    server.serve_forever()

def main(global_config, **settings):
    import os
    engine = engine_from_config(settings, 'sims.')
    listen(engine, "connect", do_connect)
    listen(engine, "begin", do_begin)
    DBSession.configure(bind=engine)
    Base.metadata.bind = engine
    Base.metadata.create_all(engine)
    config = Configurator(settings=settings)
    config.scan()
    config.include('pyramid_chameleon')

    here = os.path.dirname(os.path.abspath(__file__))
    config.add_route('huelsenbeck_svg', '/huelsenbeck_svg/{tree_estimate_set_id}.svg')
    config.add_view(huelsenbeck_svg, route_name='huelsenbeck_svg')
    #config.add_view(huelsenbeck_svg, route_name='huelsenbeck_svg', renderer=os.path.join(here,'templates','huelsenbeck_tree_estimate_set.pt'))
    config.add_route('huelsenbeck_tree_estimate_set', '/huelsenbeck_tree_estimate_set/{tree_estimate_set_id}')
    config.add_view(huelsenbeck_tree_estimate_set, route_name='huelsenbeck_tree_estimate_set', renderer=os.path.join(here,'templates','huelsenbeck_tree_estimate_set.pt'))
    config.add_route('huelsenbeck', '/huelsenbeck')
    config.add_route('huelsenbeck_tree_estimate_sets', '/huelsenbeck_tree_estimate_sets/{method}/{distance_formula}/{k}')
    config.add_route('huelsenbeck_simulation_sets', '/huelsenbeck_simulation_sets/{n}/{indelible_model}/{genes}/{m}/{theta}')
    config.add_view(huelsenbeck, route_name='huelsenbeck', renderer=os.path.join(here,'templates','huelsenbeck.pt'))
    config.add_view(huelsenbeck, route_name='huelsenbeck_tree_estimate_sets', renderer=os.path.join(here,'templates','huelsenbeck.pt'))
    config.add_view(huelsenbeck, route_name='huelsenbeck_simulation_sets', renderer=os.path.join(here,'templates','huelsenbeck.pt'))
    config.add_route('results', '/results')
    config.add_view(results, route_name='results', renderer=os.path.join(here,'templates','results.pt'))
    static_dir = settings['static_dir']
    config.add_static_view(name='static', path='kmers:static/')
    config.override_asset(to_override='kmers:static/',
                                  override_with=static_dir)
    return config.make_wsgi_app()
