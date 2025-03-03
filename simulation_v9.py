#!/usr/bin/env python3
# coding=utf-8
import sys
import time
import argparse
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import collections as mc
from matplotlib import cm
import copy
import pdb
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
from matplotlib import colors
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
matplotlib.rcParams.update({'font.size': 6})

def get_parser():
    # build commandline parser
    parser = argparse.ArgumentParser(
        description="Script for simulating hyphal growth.")

    # arguments currently not in use - parameters are hardcoded in the script
    # arguments could easily be integrated for easier tests of parameter changes (in bash)
    # I have included all parameters that we can adjust
    """
    parser.add_argument("-minTheta", type=int, dest="minTheta", default=10,
                        help="*pi/180 - minimum angle a branch can occur. Default: 10")
    parser.add_argument("-maxTheta", type=int, dest="maxTheta", default=100,
                        help="*pi/180 - maximum angle a branch can occur. Default: 100")
    parser.add_argument("-avgRate", type=int, dest="avgRate", default=50,
                        help="average tip extension rate (?m/h), source: Spohr et al. 1998. Default: 50")
    parser.add_argument("-q", type=float, dest="q", default=0.005,
                        help="branching frequency (maybe need to be scaled of what the time step is). Default: 0.005")
    parser.add_argument("-N", type=int, dest="N", default=1000,
                        help="max simulation rounds. Default: 1000")
    parser.add_argument("-M", type=int, dest="M", default=100000,
                        help="Max hyphae (length?). Default: 100000")
    parser.add_argument("-tstep", type=float, dest="tstep", default=0.0005,
                        help="Time step (h). Default: 0.0005")
    parser.add_argument("-ktip1", type=int, dest="ktip1", default=80,
                        help="Initial tip extension rate, value estimated from Spohr et al 1998 figure 5. Default: 80")
    parser.add_argument("-maxktip", type=int, dest="maxktip", default=155,
                        help="Maximum tip extension rate, value estimated from Spohr et al 1998 figure 5. Default: ")
    parser.add_argument("-Kt", type=int, dest="Kt", default=5,
                        help="Saturation constant. Default: 5")
    """
    return parser

def get_args():
    parser = get_parser()
    args = parser.parse_args()
    return args

def timepoint(comment, timelist):
    now = time.time()
    timelist.append((comment, now))

def showtimepoints(timelist):
    m = max([ len(timelist[i][0]) for i in range(1, len(timelist)) ])
    formatstring = "{:" + str(m+1) +"} {:.02f} seconds"
    for i in range(1, len(timelist)):
        print(formatstring.format(timelist[i][0]+':', timelist[i][1]-timelist[i-1][1]))
    print(formatstring.format('Total:', timelist[-1][1]-timelist[0][1]))

def plot_initial_hyphae(hyphae, tMax=100, xlim=(-100,100), ylim=(-100,100)):
    """In R, you get the points if coordinates connecting line segments are identical when doing normal scatterplot.
    But with pyplot you don't - therefore this function for plotting initial hyphae"""
    fig, ax = plt.subplots(figsize=(5, 5))
    points = list(zip(hyphae['x0'],hyphae['y0']))
    plt.plot(points, ',', color='black')
    plt.show()

def plot_hyphae(hyphae, tMax=100, xlim=(-100,100), ylim=(-100,100)):
    """Works as function in R-script, but currently not in use"""
    fig, ax = plt.subplots(figsize=(5, 5))
    line_segments = [[list(zip(hyphae['x0'],hyphae['y0']))[i], list(zip(hyphae['x'],hyphae['y']))[i]] for i in range(len(hyphae))]
    lc = mc.LineCollection(line_segments, linewidths=2)
    ax.add_collection(lc)
    ax.autoscale()
    ax.margins(0.1)
    plt.show()

def grid_bounding_boxes(w=10, xrng=(-70, 70), yrng=(-70, 70)):
    """Function making the same grid_bounding_boxes as in R-script"""
    x0 = [no for no in range(xrng[0], xrng[1]+w, w)]
    y0 = [no for no in range(yrng[0], yrng[1]+w, w)]

    xcenters, ycenters = list(), list()
    for i in range(len(x0)):
        xcenters.extend(x0)
        ycenters.extend([y0[i]]*len(y0))
    centers = pd.DataFrame({0:xcenters, 1:ycenters})
    bbs = pd.concat([centers-(w/2), centers.rename(columns={0:2, 1:3})+(w/2)], axis=1)

    return bbs

def grid_bounding_boxes_lower_left(w=10, xrng=(-70, 70), yrng=(-70, 70)):
    """grid_bounding_boxes function for plotting substrate rectangles
    Returns a dataframe of coordinates of all lower left corners"""
    x0 = [no for no in range(xrng[0], xrng[1]+w, w)]
    y0 = [no for no in range(yrng[0], yrng[1]+w, w)]

    xcenters, ycenters = list(), list()
    for i in range(len(x0)):
        xcenters.extend(x0)
        ycenters.extend([y0[i]]*len(y0))
    centers = pd.DataFrame({0:xcenters, 1:ycenters})

    return centers-(w/2)

def hyphae_hits_substrate(hl, bbs, hyph_idx):
    m = len(hl)
    h2b = pd.DataFrame(np.zeros((m, len(bbs))))

    for hl_row in range(m):
        x_val = hl[hl_row][hyph_idx['x']]
        y_val = hl[hl_row][hyph_idx['y']]

        xSAT = (x_val <= bbs[2]) & (bbs[0] <= x_val)
        ySAT = (y_val <= bbs[3]) & (bbs[1] <= y_val)
        bbInds = (xSAT & ySAT).index[(xSAT & ySAT) == True]   # index of the True in boolean series
        h2b.iloc[hl_row, bbInds] = 1
    return h2b

def tipExtensionMonod(ktip1, ktip2, Kt, l, S, Ks):
    """source: Lejeune et al 1995, Morphology of Trichoderma reesei QM 9414 in Submerged Cultures"""
    return (ktip1 + ktip2*(l/(l+Kt))) * S/(S+Ks)

def hyphal_length(x, hyph_idx):
  return math.sqrt( (x[hyph_idx['x']]-x[hyph_idx['x0']])**2 + (x[hyph_idx['y']]-x[hyph_idx['y0']])**2 )

def main(args):
    # PARAMETERS USED IN THE SIMULATION
    leftright = [-1, 1]  # branching left or right
    minTheta = 10  # *pi/180 - minimum angle a branch can occur
    maxTheta = 100  # *pi/180 - maximum angle a branch can occur

    avgRate = 50  # average tip extension rate (?m/h), source: Spohr et al. 1998
    q = 0.005  # branching frequency (maybe need to be scaled of what the time step is)
    N = 200  # max simulation rounds
    M = 100000  # max hyphae
    tstep = 0.0005  # time step (h)

    ktip1 = 80  # initial tip extension rate, value estimated from Spohr et al 1998 figure 5
    maxktip = 155  # maximum tip extension rate, value estimated from Spohr et al 1998 figure 5
    ktip2 = maxktip - ktip1  # difference between ktip1 and maxktip
    Kt = 5  # saturation constant

    # Single square as growth area
    # put init_n spores on centers and define maximum width
    centers = [[0,0]]
    width = 50
    init_n = 10

    # lists for saving results during the simulation
    snapShots = dict()
    hyph_idx = {'x0':0,'y0':1, 'x':2,'y':3,'angle':4, 'nevents':5, 't':6,'l':7}
    #hyphae_dict = {'x0':[],'y0':[], 'x':[],'y':[],'angle':[]}
    hyphae = []

    hyphaeSnapshots = dict()

    # initializes the first spore(s) and growth angles
    hNum = 0
    for i in range(len(centers)):
        for j in range(init_n):
            rxy = np.random.uniform(-width, width, 2) + [int(no) for no in centers[i]]
            iTheta = int(np.random.uniform(0,360))

            # initialize spores
            hyphae.append([rxy[0],rxy[1],rxy[0]+0.01,rxy[1]+0.01,iTheta,0,0,0])
            hNum += 1

    # append initial dataframe to snapshot lists
    snapShots[0] = pd.DataFrame(copy.deepcopy(hyphae), columns=('x0','y0','x','y','angle','nevents','t','l'))
    hyphaeSnapshots['0'] = pd.DataFrame(copy.deepcopy(hyphae),columns=('x0','y0','x','y','angle','nevents','t','l'))
    #plot_initial_hyphae(snapShots[0])

    S0 = 50000                      # initial substrate concentration
    St = [S0 for i in range(N)]     # track total substrate level over time
    S_snapshots = dict()            # substrate snapshots

    # substrate grid setup
    S_bbs = grid_bounding_boxes(w=10, xrng=(-70, 70), yrng=(-70, 70))
    S = pd.DataFrame({'S':[S0/len(S_bbs) for i in range(len(S_bbs))]})
    S_snapshots[0] = S

    dls = [0]*N                     # track biomass densities over time
    r = 1                           # forgot specific name, but used to scale biomass density

    timepoint('Initializing spores and data handlers', timelist)

    # actual simulation
    i, m = 0, 0
    while i < N and m < M:
        m = len(hyphae)             # number of hyphae at this step
        #pdb.set_trace()
        # get where each hyphae tip is in the substrate grid
        hyphae_substrate = hyphae_hits_substrate(hl=hyphae, bbs=S_bbs, hyph_idx=hyph_idx)      # takes some time
        #timepoint('hyphae_hits_substrate', timelist)

        dl = 0                      # track biomass at this step

        # loop through hyphae, find new coordinates for tips
        for index in range(m):
            hyphae_row=hyphae[index]

            # check if enough substrate around substrate tip (>0)
            try:
                if hyphae_substrate.iloc[index,:].sum():
                    angle = hyphae_row[hyph_idx['angle']]
                    hi = [hyphae_row[hyph_idx['x0']], hyphae_row[hyph_idx['y0']], hyphae_row[hyph_idx['x']], hyphae_row[hyph_idx['y']]]     # OBS! no need to include x0 and y0 - not used (neither in R code)

                    # tip extension value, calculated as the accelerated growth * substrate limitation
                    # equation of the Monod type
                    extension = tipExtensionMonod(ktip1=ktip1, ktip2=ktip2, Kt=Kt, l=hyphae_row[hyph_idx['l']],
                                                  S=S.loc[hyphae_substrate.iloc[index,:] == 1].iloc[0,0], Ks=200)
                    dx = extension * tstep * math.cos(angle * math.pi / 180)    # new coordinate in x-axis
                    dy = extension * tstep * math.sin(angle * math.pi / 180)    # new coordinate in y-axis
                    #timepoint('Tip extension Monod', timelist)

                    # biomass just created for hyphae j
                    dl_c = math.sqrt(dx**2 + dy**2)
                    dl += dl_c

                    # substrate used by the hyphal tip
                    S.loc[hyphae_substrate.iloc[index, :] == 1] = S.loc[hyphae_substrate.iloc[index, :] == 1].iloc[0, 0] - min(r*dl_c, S.loc[hyphae_substrate.iloc[index, :] == 1].iloc[0, 0])
                    #timepoint('substrate used by the hyphal tip', timelist)

                    # update data frame
                    hyphae[index][hyph_idx['x']] = hi[2] + dx
                    hyphae[index][hyph_idx['y']] = hi[3] + dy
                    hyphae[index][hyph_idx['l']] = hyphal_length(hyphae_row, hyph_idx)
                    #b = hyphae['nevents'].iloc[j]
                    #timepoint('Update df', timelist)

                    qApl = q
                    # qApl = q/(b+1)        # q could be a function of number of events
                    # branching event (q should depend on l and/or nevents)

                    if np.random.uniform(0,1) < qApl:
                        # the direction a new branch will grow - randomly left or right
                        newdir = np.random.choice(leftright)
                        newangle = angle + newdir * round(np.random.uniform(minTheta,maxTheta))

                        # nevents goes one up, as there have been branching event for this hyphae
                        hyphae[index][hyph_idx['nevents']] += 1

                        # add new hyphae
                        new_hyphae = [hi[2], hi[3], hi[2], hi[3], newangle, 0, i, 0]
                        hyphae.append(new_hyphae)
            except:
                pdb.set_trace()

        # subtract the resources consumed in the last round of growth
        # S[i+1] = S[i] - min(r*dl, S[i])
        St[i] = S.sum()
        # save biomass at this timestep
        dls[i] = dl

        # increment i before saving time information
        i += 1
        if (i) % 10 == 0:
            #pdb.set_trace()
            print(f"Iteration = {i}, hNum = {len(hyphae)}")
            snapShots[i] = pd.DataFrame(copy.deepcopy(hyphae), columns=('x0','y0','x','y','angle','nevents','t','l'))
            S_snapshots[i] = copy.deepcopy(S)
            hyphaeSnapshots[str(i)] = pd.DataFrame(copy.deepcopy(hyphae), columns=('x0','y0','x','y','angle','nevents','t','l'))

    timepoint('Actual simulation', timelist)
    print('# Simulation done! Plotting...')

    """
    # Outcommented in the R code
    m = len(hyphae)
    dat = copy.deepcopy(hyphae)
    x_range = (dat[['x0', 'x']].min().min(), dat[['x0', 'x']].max().max())
    y_range = (dat[['y0', 'y']].min().min(), dat[['y0', 'y']].max().max())
    tMax = max(dat['t'])
    """

    # dimensions for subplots on one page (n-rows and m-cols)
    n, m = 3, 4
    pdf_name = 'output_test_test.pdf'
    with PdfPages(pdf_name) as pdf:
        # initialize layout for plots
        f, axarr = plt.subplots(n, m, sharex='none', sharey='none')
        arr_ij = [(x, y) for x, y in np.ndindex(axarr.shape)]
        subplots = [axarr[index] for index in arr_ij]

        splot_index = 0

        # iterate through snapshots, selecting every second
        for i in range(1, len(snapShots), 2):
            # select every second snapshot sorted by time
            snapshot_to_plot = snapShots[sorted(snapShots)[i]]
            S_snapshot_to_plot = S_snapshots[sorted(snapShots)[i]]
            line_segments = [[list(zip(snapshot_to_plot['x0'], snapshot_to_plot['y0']))[i], list(zip(snapshot_to_plot['x'], snapshot_to_plot['y']))[i]] for i in
                             range(len(snapshot_to_plot))]

            # create continuous norm for mapping colors to hyphae according to time they occured
            hyphae_norm = plt.Normalize(min(hyphae[hyph_idx['t']]), max(hyphae[hyph_idx['t']]))
            lc = mc.LineCollection(line_segments, linewidths=0.5, cmap='GnBu_r', norm=hyphae_norm)
            lc.set_array(snapshot_to_plot['t'])
            subplots[splot_index].add_collection(lc)

            # create continuous norm for substrate color mapping
            substrate_norm = plt.Normalize(0, S0/len(S_bbs))
            #subplots[splot_index].imshow(np.array(S).reshape(15,15), cmap='YlOrBr_r', norm=substrate_norm)
            S_plot_bbs = grid_bounding_boxes_lower_left(w=10, xrng=(-70, 70), yrng=(-70, 70))

            YlOrBr = cm.get_cmap('YlOrBr')
            patches, colors = list(), list()
            for j in range(len(S_plot_bbs)):
                rect = plt.Rectangle((S_plot_bbs.iloc[j,0], S_plot_bbs.iloc[j,1]), width=10, height=10) #color=YlOrBr(S_snapshot_to_plot.iloc[j]/S0/len(S_bbs))
                patches.append(rect)
            #pdb.set_trace()

            p = PatchCollection(patches, cmap=cm.get_cmap('YlOrBr'), alpha=0.5, norm=substrate_norm)
            p.set_array(S_snapshot_to_plot['S'])
            subplots[splot_index].add_collection(p)

            subplots[splot_index].autoscale()
            subplots[splot_index].margins(0.1)
            subplots[splot_index].set_title(f"Time = {round(sorted(snapShots)[i]*tstep,3)}")
            subplots[splot_index].set_xlim(-60, 60)
            subplots[splot_index].set_ylim(-60, 60)
            subplots[splot_index].set_xticks(list(range(-60,70,30)))
            subplots[splot_index].set_yticks(list(range(-60,70,30)))
            plt.tight_layout()

            # increment index counter
            splot_index += 1
            # a pdf page is full once the count is equal to n*m
            if splot_index == n * m:
                pdf.savefig()
                plt.close(f)
                f, axarr = plt.subplots(n, m, sharex='none', sharey='none')
                arr_ij = [(x, y) for x, y in np.ndindex(axarr.shape)]
                subplots = [axarr[index] for index in arr_ij]
                splot_index = 0

        # save last page
        pdf.savefig()
        plt.close(f)
        timepoint('Plotting to pdf', timelist)
    return pdf_name


if __name__ == '__main__':
    timelist = list()
    timepoint('Start', timelist)
    print("# Running...")
    args = get_args()
    print("# args:", args)
    pdf_name = main(args)
    print(f'# Done! Made pdf: {pdf_name}')
    # Display timing
    showtimepoints(timelist)
