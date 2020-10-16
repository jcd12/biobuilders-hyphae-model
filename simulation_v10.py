#!/usr/bin/env python3
# coding=utf-8
import sys
import time
import argparse
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import collections as mc
import copy
import pdb
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
from matplotlib import colors
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import imageio as io
import os
import matplotlib.colors as clr

matplotlib.rcParams.update({'font.size': 6})

def get_parser():
    # build commandline parser
    parser = argparse.ArgumentParser(
        description="Script for simulating hyphal growth.")
    parser.add_argument("-img", "--img_ana", type=str, dest="img_ana",
                        help="Output csv from image analysis tool")
    parser.add_argument("-mu_max", type=float, dest="mu_max", default=0.275733333333333,
                        help="Mu_max. Maximum growth rate for the given strain.")

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

def img_ana_parameters(csv_file):
    return pd.read_csv(csv_file, header=0)

def plot_initial_hyphal_elements(hyphal_elements, tMax=100, xlim=(-100,100), ylim=(-100,100)):
    """In R, you get the points if coordinates connecting line segments are identical when doing normal scatterplot.
    But with pyplot you don't - therefore this function for plotting initial hyphal_elements"""
    fig, ax = plt.subplots(figsize=(5, 5))
    points = list(zip(hyphal_elements['x0'],hyphal_elements['y0']))
    plt.plot(points, ',', color='black')
    plt.show()

def plot_hyphal_elements(hyphal_elements, tMax=100, xlim=(-100,100), ylim=(-100,100)):
    """Works as function in R-script, but currently not in use"""
    fig, ax = plt.subplots(figsize=(5, 5))
    line_segments = [[list(zip(hyphal_elements['x0'],hyphal_elements['y0']))[i], list(zip(hyphal_elements['x'],hyphal_elements['y']))[i]] for i in range(len(hyphal_elements))]
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

def tip_extension_monod(mu_max, sep_row, hyphal_elements, St, Ks, extension_substrate_min, curvature_gamma_params, r, time):
    """source: Lejeune et al 1995, Morphology of Trichoderma reesei QM 9414 in Submerged Cultures"""
    sep_distance = center_distance(sep_row)
    if sep_distance >= 20:
        return hyphal_elements

    extension = mu_max * St[int(r/20*sep_distance)] / (St[int(r/20*sep_distance)] + Ks)
    if St[int(r/20*sep_distance)] > extension_substrate_min and extension > np.random.uniform(0, 1):
        # get coordinates of old tip
        x0_old, y0_old, x1_old, y1_old = midpoints_to_endpoints(sep_row)
        if 0 < sep_row.angle < 180:
            x_oldtip, y_oldtip = max(((x0_old, y0_old), (x1_old, y1_old)),  key=lambda x: x[1])
        elif 180 < sep_row.angle < 360:
            x_oldtip, y_oldtip = min(((x0_old, y0_old), (x1_old, y1_old)),  key=lambda x: x[1])
        elif sep_row.angle == 0:
            x_oldtip, y_oldtip = max(((x0_old, y0_old), (x1_old, y1_old)), key=lambda x: x[0])
        elif sep_row.angle == 180:
            x_oldtip, y_oldtip = min(((x0_old, y0_old), (x1_old, y1_old)), key=lambda x: x[0])

        new_coord, center_dist, angles = list(), list(), list()
        for _ in range(10):
            new_angle = sep_row.angle + np.random.choice((-1, 1)) * round(5*np.random.gamma(*curvature_gamma_params))
            # check that new angle is >=0 and < 360
            if not 0 <= new_angle < 360:
                new_angle = new_angle % 360

            # find out if x and y positions for new branch are smaller/larger than before
            if 0 < new_angle < 180:
                y_dir = 1
            else:
                y_dir = -1
            if 0 < new_angle < 90 or 270 < new_angle < 360:
                x_dir = 1
            else:
                x_dir = -1

            # make sure you get right angle of triangle
            if 0 < new_angle < 90: use_angle = new_angle
            elif 90 < new_angle < 180:
                use_angle = 180 - new_angle
            elif 180 < new_angle < 270:
                use_angle = new_angle - 180
            elif 270 < new_angle < 360:
                use_angle = 360 - new_angle

            if new_angle in (0, 90, 180, 270):
                if new_angle in (0, 180):
                    delta_x = x_dir * 0.5
                    delta_y = 0
                elif new_angle in (90, 270):
                    delta_x = 0
                    delta_y = y_dir * 0.5

            else:
                # change in y direction based on sinus relations
                delta_y = y_dir * abs(np.sin(np.deg2rad(use_angle)) / (2 * np.sin(np.deg2rad(90))))
                # change in x direction based on normal pythagoras
                delta_x = x_dir * np.sqrt(0.5 ** 2 - delta_y ** 2)

            x_mid, y_mid = x_oldtip + delta_x, y_oldtip + delta_y

            new_coord.append((x_mid, y_mid))
            dist = np.sqrt(x_mid ** 2 + y_mid ** 2)
            if dist > 20: dist = 20
            center_dist.append(dist)
            angles.append(new_angle)

            # if substrate concentration is more than 1% as high for one proposed angle than the first,
            # choose this angle, else choose first
            if St[int(r/20*max(center_dist))] > 1.01 * St[int(r/20*center_dist[0])]:
                print('tropism!')
                new_x_mid, new_y_mid = new_coord[np.argmax(center_dist)]
                new_angle = angles[np.argmax(center_dist)] % 360

            else:
                new_x_mid, new_y_mid = new_coord[0]
                new_angle = angles[0] % 360

        
        # OBS!: make new angle dependent on substrate concentration - atm 10 angles are proposed and the one maximising distance to center is chosen
        """
        new_angle = sep_row.angle + np.random.choice((-1,1)) * round(np.random.uniform(min_curvature, max_curvature))
        # check that new angle is >=0 and < 360
        if not 0 <= new_angle < 360:
            new_angle = new_angle % 360

        # find out if x and y positions for new branch are smaller/larger than before
        if 0 < new_angle < 180: y_dir = 1
        else: y_dir = -1
        if 0 < new_angle < 90 or 270 < new_angle < 360: x_dir = 1
        else: x_dir = -1

        # make sure you get right angle of triangle
        if 0 < new_angle < 90: use_angle = new_angle
        if 90 < new_angle < 180: use_angle = 180 - new_angle
        elif 180 < new_angle < 270: use_angle = new_angle - 180
        elif 270 < new_angle < 360: use_angle = 360 - new_angle

        if new_angle in (0, 90, 180, 270):
            if new_angle in (0,180):
                delta_x = x_dir*0.5
                delta_y = 0
            elif new_angle in (90,270):
                delta_x = 0
                delta_y = y_dir*0.5

        else:
            # change in y direction based on sinus relations
            delta_y = y_dir * abs(np.sin(np.deg2rad(use_angle)) / (2*np.sin(np.deg2rad(90))))
            # change in x direction based on normal pythagoras
            delta_x = x_dir * np.sqrt(0.5**2-delta_y**2)

        new_x_mid, new_y_mid = x_oldtip + delta_x, y_oldtip + delta_y


        
        proposed_angles = [sep_row.angle + np.random.choice((-1, 1)) * round(np.random.uniform(min_curvature, max_curvature)) for _ in range(10)]
        center_dist, new_coord = list(), list()
        for angle in proposed_angles:
            # check that new angle is >=0 and < 360
            if not 0 <= angle < 360: angle = angle % 360

            # find out if x and y positions for new branch are smaller/larger than before
            if 0 < angle < 180: y_dir = 1
            else: y_dir = -1
            if 0 < angle < 90 or 270 < angle < 360: x_dir = 1
            else: x_dir = -1

            # make sure you get right angle of triangle
            if 0 < angle < 90: use_angle = angle
            if 90 < angle < 180: use_angle = 180 - angle
            elif 180 < angle < 270: use_angle = angle - 180
            elif 270 < angle < 360: use_angle = 360 - angle

            if angle in (0, 90, 180, 270):
                if angle in (0, 180):
                    delta_x = x_dir * 0.5
                    delta_y = 0
                elif angle in (90, 270):
                    delta_x = 0
                    delta_y = y_dir * 0.5
            else:
                # change in y direction based on sinus relations
                delta_y = y_dir * abs(np.sin(np.deg2rad(use_angle)) / (2 * np.sin(np.deg2rad(90))))
                # change in x direction based on normal pythagoras
                delta_x = x_dir * np.sqrt(0.5 ** 2 - delta_y ** 2)

                new_coord.append((x_oldtip + delta_x, y_oldtip + delta_y))

                # distance to center
                center_dist.append(np.sqrt((x_oldtip + delta_x) ** 2 + (y_oldtip + delta_y) ** 2))

        # if substrate concentration is more than twice as high for one proposed angle than the first,
        # choose this angle, else choose first
        if St[round(max(center_dist))] > 1.5 * St[round(center_dist[0])]:
            print('tropism!')
            new_x_mid, new_y_mid = new_coord[np.argmax(center_dist)]
            new_angle = proposed_angles[np.argmax(center_dist)] % 360
        else:
            new_x_mid, new_y_mid = new_coord[0]
            new_angle = proposed_angles[0] % 360
        """
        hyphal_elements.loc[len(hyphal_elements)] = {'x_mid':new_x_mid, 'y_mid':new_y_mid, 'angle': new_angle, 'tip': True, 'time':time}

        # change tip status of sep_row, as this hyphal_elements has been extended and is no longer a tip
        hyphal_elements.at[sep_row.name, 'tip'] = False
    return hyphal_elements

def hyphal_elements_hits_substrate(hl, bbs):
    """
    :param hl:
    :param bbs:
    :return: dataframe with row for each hyphal_elements of length number of
    bounding boxes - only one 1 per row telling you which bounding
    box the given tip is in
    """
    m = len(hl)
    h2b = pd.DataFrame(np.zeros((m, len(bbs))))

    def get_hyphal_elements_position(row):
        idx = row.name
        x_val = row['x_mid']
        y_val = row['y_mid']

        xSAT = (x_val <= bbs[2]) & (bbs[0] <= x_val)
        ySAT = (y_val <= bbs[3]) & (bbs[1] <= y_val)
        bbInds = (xSAT & ySAT).index[(xSAT & ySAT) == True]     # index of the True in boolean series
        h2b.iloc[idx, bbInds] = 1
        return

    hl.apply(lambda row: get_hyphal_elements_position(row), axis=1)

    return h2b

def hyphal_length(hyphal_elements):
  return math.sqrt( (hyphal_elements.x-hyphal_elements.x0)**2 + (hyphal_elements.y-hyphal_elements.y0)**2)

def substrate_per_distance(n_tip, n_nontip, St, S0=200, r=200, h=1, D=0.1, S_tip=5, S_nontip=5):
    """
    :param n_tip: number of tip hyphal_elements
    :param n_nontip: number of non-tip hyphal_elements
    :param St: St list of substrate values for previous simulation step
    :param S0: initial substrate concentration - also boundary condition
    :param r: radius of substrate field
    :param h: radial step size
    :param D: diffusion constant
    :param S_tip: amount of substrate used per tip septum
    :param S_nontip: amount of substrate used per non-tip septum
    :return: St: list of increasing substrate concentrations as a function of distance to center
    """

    dSdt = [0 for i in range(r + 1)]  # list to hold changes in substrate levels for this timestep
    # continuous equation is dS/dt = D*d2S/dr2 (minus use at centre) - using difference estimates for second order derivative
    # forward difference for center
    dSdt[0] = min(0, D * (St[0 + 2 * h] - 2 * St[1] + St[0]) / h ** 2 - (S_tip * n_tip + S_nontip * n_nontip))
    for rad in range(1, r):
        # centre difference for most of the field - also discarding positive differences which happen due to numerical solution
        dSdt[rad] = min(0, D * (St[rad + 1] - 2 * St[rad] + St[rad - 1]) / h ** 2)
    # centre difference but with S(rad>r)=S_naught
    dSdt[r] = min(0, D * (-St[r] + St[r - 1]) / h ** 2)
    # OBS!: change to 

    #St[0] = St[0] + dSdt[0]
    St = [max(0, a + b) for a, b in zip(St, dSdt)]

    # raise error since I don't think local concentrations above S0 makes sense
    # OBS!: this is commonly detected - printed to test

    if max(St) > S0:
        St = [min(st, S0) for st in St]
        print("instability detected")

    # OBS!: something wrong with first pos - often large
    #St[0] = 0
    return St


def center_distance(sep_row):
    return np.sqrt(sep_row.x_mid**2+sep_row.y_mid**2)

def branching(sep_row, St, min_substrate, p_branching, hyphal_elements, branch_angle_beta_params, r, time):
    """
    Add laterally or apically branched hyphal_elements (by a certain probability) to hyphal_elements df if enough substrate
    :param sep_row: row from hyphal_elements df
    :param St: substrate concentration per distance for this timestep
    :param min_substrate: minimum substrate for branching to occur
    :param p_branching: probability of branching event, given enough substrate
    :param hyphal_elements: hyphal_elements df to add evt. branched septum to
    :return:
    """
    # OBS!: right now apical branching also occurs from the middle of the given hyphal_elements!
    # if this should be changed, we should consider whether given branching event is
    # apical or lateral and then for apical, first calculate the position of the tip
    # but this requires information on which end of the hyphal_elements is the tip-end (can be seen from the angle?)
    sep_distance = center_distance(sep_row)
    if sep_distance >= 20:
        return hyphal_elements

    # lateral branching occurs at a given probability if enough substrate
    if St[int(r/20*sep_distance)] > min_substrate and p_branching > np.random.uniform(0, 1):

        print('BRANCHING!')
        new_angle = sep_row.angle + np.random.choice((-1,1)) * round(np.random.beta(*branch_angle_beta_params)*90)      # 90 is the largest branching angle
        # check that new angle does not exceed 360
        if not 0 <= new_angle < 360:
            new_angle = new_angle % 360

        # find out if x and y positions for new branch are smaller/larger than before
        if 0 < new_angle < 180: y_dir = 1
        else: y_dir = -1
        if 0 < new_angle < 90 or 270 < new_angle < 360: x_dir = 1
        else: x_dir = -1

        # make sure you get right angle of triangle
        if 0 < new_angle < 90: use_angle = new_angle
        if 90 < new_angle < 180: use_angle = 180 - new_angle
        elif 180 < new_angle < 270: use_angle = new_angle - 180
        elif 270 < new_angle < 360: use_angle = 360 - new_angle

        if new_angle in (0, 90, 180, 270):
            if new_angle in (0,180):
                delta_x = x_dir*0.5
                delta_y = 0
            elif new_angle in (90,270):
                delta_x = 0
                delta_y = y_dir*0.5

        else:
            # change in y direction based on sinus relations
            delta_y = y_dir * abs(np.sin(np.deg2rad(use_angle)) / (2*np.sin(np.deg2rad(90))))
            # change in x direction based on normal pythagoras
            delta_x = x_dir * np.sqrt(0.5**2-delta_y**2)

        new_x_mid, new_y_mid = sep_row.x_mid + delta_x, sep_row.y_mid + delta_y

        hyphal_elements.loc[len(hyphal_elements)] = {'x_mid':new_x_mid, 'y_mid':new_y_mid, 'angle': new_angle, 'tip': True, 'time':time}

    return hyphal_elements


def midpoints_to_endpoints(sep_row):
    """Use x_mid, y_mid and angle to find end coordinates of hyphal_elements"""

    # return for easy angles
    if sep_row.angle in (0, 90, 180, 270):
        if sep_row.angle in (0,180):
            return sep_row.x_mid-0.5, sep_row.y_mid, sep_row.x_mid+0.5, sep_row.y_mid
        elif sep_row.angle in (90,270):
            return sep_row.x_mid, sep_row.y_mid - 0.5, sep_row.x_mid, sep_row.y_mid + 0.5

    # get right angle of triangle to use for sinus relations
    if sep_row.angle < 90: use_angle = sep_row.angle
    if 90 < sep_row.angle < 180: use_angle = 180 - sep_row.angle
    if 180 < sep_row.angle < 270: use_angle = sep_row.angle - 180
    if 270 < sep_row.angle < 360: use_angle = 360 - sep_row.angle

    # change in y direction based on sinus relations
    delta_y = abs(np.sin(np.deg2rad(use_angle)) / (2 * np.sin(np.deg2rad(90))))
    # change in x direction based on normal pythagoras
    delta_x = np.sqrt(abs(0.5 ** 2 - delta_y ** 2))

    # check how to add delta_x and delta_y
    if 0 < sep_row.angle < 90 or 180 < sep_row.angle < 270:
        x0, y0 = sep_row.x_mid - delta_x, sep_row.y_mid - delta_y
        x1, y1 = sep_row.x_mid + delta_x, sep_row.y_mid + delta_y
    else:
        x0, y0 = sep_row.x_mid - delta_x, sep_row.y_mid + delta_y
        x1, y1 = sep_row.x_mid + delta_x, sep_row.y_mid - delta_y

    return x0, y0, x1, y1


def plot_to_pdf(snapshots, St_snapshots, S0, r, tstep):
    # dimensions for subplots on one page (n-rows and m-cols)
    n, m = 3, 4
    pdf_name = 'ATCC_12h.pdf'
    with PdfPages(pdf_name) as pdf:
        # initialize layout for plots
        f, axarr = plt.subplots(n, m, sharex='none', sharey='none')
        arr_ij = [(x, y) for x, y in np.ndindex(axarr.shape)]
        subplots = [axarr[index] for index in arr_ij]

        splot_index = 0

        S_max_index = np.argmax(St_snapshots[max(St_snapshots)])

        for i in range(0, len(snapshots)):
            snapshot_to_plot = snapshots[sorted(snapshots)[i]]
            St_to_plot = St_snapshots[sorted(snapshots)[i]]

            endpoints_df = pd.DataFrame({'x0y0x1y1': snapshot_to_plot.apply(midpoints_to_endpoints, axis=1)})
            x0, y0, x1, y1 = endpoints_df['x0y0x1y1'].str[0].to_list(), endpoints_df['x0y0x1y1'].str[1].to_list(), \
                             endpoints_df['x0y0x1y1'].str[2].to_list(), endpoints_df['x0y0x1y1'].str[3].to_list()

            line_segments = [[list(zip(x0, y0))[i], list(zip(x1,y1))[i]] for i in
                             range(len(snapshot_to_plot))]
            # create continuous norm for mapping colors to hyphae according to time they occured
            hyphal_elements_norm = plt.Normalize(-10, len(snapshots)*1.25)
            red_cmap = clr.LinearSegmentedColormap.from_list('custom red', [(0, '#690909'), (0.5, '#AE1212'), (0.75, '#FD8D4C'),
                                                                (1, '#FCE0CB')],N=400)
            lc = mc.LineCollection(line_segments, linewidths=0.5, norm=hyphal_elements_norm, cmap=red_cmap) # OrRd_r
            lc.set_array(snapshot_to_plot['time'])
            subplots[splot_index].add_collection(lc)

            substrate_norm = plt.Normalize(0, S0)
            xlist = np.linspace(-20, 20, 100)
            ylist = np.linspace(-20, 20, 100)
            X, Y = np.meshgrid(xlist, ylist)
            dist = np.sqrt(X ** 2 + Y ** 2)
            dist[dist > 20] = 20
            Z = [[St_to_plot[int(r/20*x)] for x in y] for y in dist]
            subplots[splot_index].contourf(X, Y, Z, cmap=cm.get_cmap('YlGn'), norm=substrate_norm, alpha=0.7)

            subplots[splot_index].autoscale()
            subplots[splot_index].margins(0.1)
            subplots[splot_index].set_title(
                f"Time: {round(2 * i * tstep, 2)} h\n"
                f"Number of hyphal_elements: {len(snapshot_to_plot)}\n"
                f"Number of hyphal tips: {len(snapshot_to_plot[snapshot_to_plot.tip == True])}\n"
                f"Branching frequency: {round((len(snapshot_to_plot[snapshot_to_plot.tip == True]) / len(snapshot_to_plot)) / (tstep * i + float('1e-10')), 3)}\n"
                f"(branches per hyphal element per hour)",
                fontsize=5)
            subplots[splot_index].set_xlim(-20, 20)
            subplots[splot_index].set_ylim(-20, 20)
            subplots[splot_index].axis('off')
            #subplots[splot_index].set_xticks(list(range(-60,70,30)))
            #subplots[splot_index].set_yticks(list(range(-60,70,30)))

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
    return pdf_name


def plot_for_animation(snapshots, St_snapshots, S0, dirname, r, tstep, max_time):

    for i in range(len(snapshots)):
        snapshot_to_plot = snapshots[sorted(snapshots)[i]]
        St_to_plot = St_snapshots[sorted(snapshots)[i]]

        endpoints_df = pd.DataFrame({'x0y0x1y1': snapshot_to_plot.apply(midpoints_to_endpoints, axis=1)})
        x0, y0, x1, y1 = endpoints_df['x0y0x1y1'].str[0].to_list(), endpoints_df['x0y0x1y1'].str[1].to_list(), \
                         endpoints_df['x0y0x1y1'].str[2].to_list(), endpoints_df['x0y0x1y1'].str[3].to_list()

        fig, ax = plt.subplots(figsize=(8,8))

        line_segments = [[list(zip(x0, y0))[i], list(zip(x1, y1))[i]] for i in
                         range(len(snapshot_to_plot))]

        # create continuous norm for mapping colors to hyphae according to time they occured
        hyphal_elements_norm = plt.Normalize(0, max_time)
        purple_cmap = clr.LinearSegmentedColormap.from_list('custom purple', [(0, '#440154FF'),(0.75, '#404788FF'), (0.95, '#9999FF'), (1,'#CCCCFF')], N=200)
        lc = mc.LineCollection(line_segments, linewidths=1, norm=hyphal_elements_norm,
                               cmap=purple_cmap)  # cmap='OrRd_r'
        lc.set_array(snapshot_to_plot['time'])
        ax.add_collection(lc)

        substrate_norm = plt.Normalize(0, S0)  # S0
        green_cmap = clr.LinearSegmentedColormap.from_list('custom green', ['#FFFFFF', '#55C667FF'], N=200) #73D055FF
        xlist = np.linspace(-20, 20, 1000)
        ylist = np.linspace(-20, 20, 1000)
        X, Y = np.meshgrid(xlist, ylist)
        dist = np.sqrt(X ** 2 + Y ** 2)
        dist[dist>20] = 20
        Z = [[St_to_plot[int(r/20*x)] for x in y] for y in dist]

        ax.contourf(X, Y, Z, cmap=green_cmap, norm=substrate_norm,
                                       alpha=0.5)  # norm=substrate_norm    #YlGn cm.get_cmap('YlGn')

        ax.autoscale()
        #ax.margins(0.1)
        ax.set_title(
            f"Time: {round(i*tstep,2)} h\n"
            f"Number of hyphal_elements: {len(snapshot_to_plot)}\n"
            f"Number of hyphal tips: {len(snapshot_to_plot[snapshot_to_plot.tip == True])}\n"
            f"Branching frequency: {round((len(snapshot_to_plot[snapshot_to_plot.tip == True]) / len(snapshot_to_plot))/(tstep*i+float('1e-10')),3)}\n"
            f"(branches per hyphal element per hour)",
            fontsize=9)
        ax.set_xlim(-20, 20)
        ax.set_ylim(-20, 20)
        plt.axis('off')

        plt.tight_layout()
        if not os.path.exists(dirname):
            os.mkdir(dirname)
        plt.savefig(f"{dirname}/image-{'0'*(4-len(str(i)))}{i}.png")
        plt.close(fig)

    return dirname

def make_gif(filenames, gifname):
    with io.get_writer(gifname, mode='I', duration=0.05) as writer:
        for filename in filenames:
            image = io.imread(filename)
            writer.append_data(image)
    writer.close()
    return

def main(args):
    # PARAMETERS USED IN THE SIMULATION

    # if given output from image analysis tool, use parameters from this
    if args.img_ana:
        img_ana_params = img_ana_parameters(args.img_ana)
        q = img_ana_params['branching_frequency'].iloc[0]
        curvature_gamma_params = (img_ana_params['curvature_gamma_a'].iloc[0], img_ana_params['curvature_gamma_scale'].iloc[0])
        branch_angle_beta_params = (img_ana_params['angle_beta_a'].iloc[0], img_ana_params['angle_beta_b'].iloc[0])

    # default parameters are based on image analysis of ATCC 1015
    else:
        """
        img_ana_params = ()
        q =
        curvature_gamma_params =
        branch_angle_beta_params =
        """

    #min_branch_angle, max_branch_angle = 20, 100 # *pi/180 - minimum branching angle *pi/180 - maximum branching angle
    #min_curvature = 1      # minimum curvature between extended hyphal_elements – OBS!: get from image analysis
    #max_curvature = 20     # maximum curvature between extended hyphal_elements – OBS!: get from image analysis

    avgRate = 50  # average tip extension rate (?m/h), source: Spohr et al. 1998
    q = 0.0005  # branching frequency (maybe need to be scaled of what the time step is)
    M = 2000  # max hyphal_elements

    p_lateral = 3/4*q                # OBS!: probability of lateral branching for given hyphal_elements when enough substrate – find paper!
    p_apical = 1/4*q
    lateral_substrate_min, apical_substrate_min = 10, 10      # minimum substrate level for lateral branching to occur
    extension_substrate_min = 5
    mu_max = args.mu_max            # obtained from biolector measurements

    ktip1 = 80  # initial tip extension rate, value estimated from Spohr et al 1998 figure 5
    maxktip = 155  # maximum tip extension rate, value estimated from Spohr et al 1998 figure 5
    ktip2 = maxktip - ktip1  # difference between ktip1 and maxktip
    Kt = 5  # saturation constant

    # parameters for substrate gradient
    S0 = 100    # initial substrate concentration - also boundary condition
    r = 400     # "radius" of substrate field (2000/100 = 20)
    h = 0.01       # radial step size
    St = [S0 for i in range(r + 1)]  # list of substrate values for time t as a function of distance to the centre - start with uniform field
    D = 0.99
    S_tip = 10
    S_nontip = 1
    tstep = 1/60    # 1 minute
    N = 20/tstep  # max simulation rounds
    branch_substrate_dependency = 1.3
    # lists for saving results during the simulation
    snapshots = dict()
    St_snapshots = dict()

    # initialize the first two hyphal_elements growing from a single spore
    n_hyphal_elements = 8
    hyphal_elements_dict = {'x_mid':[-0.5,0.5,0,0, np.sqrt(2)/4, -np.sqrt(2)/4, -np.sqrt(2)/4, np.sqrt(2)/4], 'y_mid':[0,0,-0.5,0.5, np.sqrt(2)/4, np.sqrt(2)/4, -np.sqrt(2)/4, -np.sqrt(2)/4], 'angle':[180,0,270,90,45,135,225,315 ], 'tip':[True, True, True, True, True, True, True, True], 'time':[0,0,0,0,0,0,0,0]}
    hyphal_elements = pd.DataFrame(hyphal_elements_dict)
    # append initial dataframe to snapshot lists
    snapshots[0] = (copy.deepcopy(hyphal_elements))
    St_snapshots[0] = St

    # actual simulation
    i, m = 1, 0
    while i < N:
        n_hyphal_elements = len(hyphal_elements)
        hyphal_elements_tip = hyphal_elements[hyphal_elements.tip == True]
        hyphal_elements_nontip = hyphal_elements[hyphal_elements.tip == False]
        print('ntip:',len(hyphal_elements_tip), 'n_nontip:', len(hyphal_elements_nontip))
        St = substrate_per_distance(n_tip=len(hyphal_elements_tip), n_nontip=len(hyphal_elements_nontip), St=St, D=D, S_tip=S_tip, S_nontip=S_nontip, r=r)
        print(St)
        # lateral branching for non-tip hyphal_elements
        for j in range(len(hyphal_elements_nontip)):
            branching(hyphal_elements_nontip.iloc[j], St, lateral_substrate_min, (branch_substrate_dependency*np.mean(St)/S0)*p_lateral, hyphal_elements, branch_angle_beta_params, r=r, time=i)

        for k in range(len(hyphal_elements_tip)):
            # apical branching for tip hyphal_elements
            branching(hyphal_elements_tip.iloc[k], St, apical_substrate_min, p_apical, hyphal_elements, branch_angle_beta_params, r=r, time=i)
            # extensions for tip hyphal_elements, using monod as extension probabilities
            tip_extension_monod(mu_max, sep_row=hyphal_elements_tip.iloc[k], hyphal_elements=hyphal_elements, St=St, Ks=200,
                                extension_substrate_min=extension_substrate_min, curvature_gamma_params=curvature_gamma_params, r=r, time=i)

        # lateral branching for non-tip hyphal_elements
        #hyphal_elements_nontip.apply(lambda sep_row: branching(sep_row, St, lateral_substrate_min, p_lateral, hyphal_elements, min_branch_angle, max_branch_angle), axis=1)
        # apical branching for tip hyphal_elements
        #hyphal_elements_tip.apply(lambda sep_row: branching(sep_row, St, apical_substrate_min, p_apical, hyphal_elements, min_branch_angle, max_branch_angle), axis=1)
        # extensions for all tip-hyphal_elementss, using monod as extension probabilities
        #hyphal_elements_tip.apply(lambda sep_row: tip_extension_monod(ktip1=ktip1, ktip2=ktip2, Kt=Kt, sep_row=sep_row, hyphal_elements=hyphal_elements, St=St, Ks=200, min_curvature=min_curvature, max_curvature=max_curvature), axis=1)

        #print(f"Iteration = {i}\tNumber of hyphal_elements = {len(hyphal_elements)}")
        #snapshots[i] = copy.deepcopy(hyphal_elements)
        #St_snapshots[i] = copy.deepcopy(St)
        i += 1
        #St_snapshots[i] = copy.deepcopy(St)
        #hyphaeSnapshots[str(i)] = copy.deepcopy(hyphae)
        
        # mangler noget? save St for plot?
        # OBS!: at some point we should check that the hyphae is not growing out of image boundaries? just check distance to center

        if (i) % 2 == 0:
            print(f"Iteration = {i}\tNumber of hyphal_elements = {len(hyphal_elements)}")
            snapshots[i] = copy.deepcopy(hyphal_elements)
            St_snapshots[i] = copy.deepcopy(St)

    pdf_name = plot_to_pdf(snapshots, St_snapshots, S0, r=r, tstep=tstep)
    #dir_name = plot_for_animation(snapshots, St_snapshots, S0, dirname='spaA_gul1_12h', r=r, tstep=tstep, max_time = i)

    #filenames = sorted(dir_name+'/'+fn for fn in os.listdir(dir_name) if fn.startswith('image'))
    #make_gif(filenames, f'{dir_name}/test.gif')
    print(snapshots[0])
    print('# Simulation done! Plotting...')
    # make function for plotting!
    return pdf_name


if __name__ == '__main__':
    timelist = list()
    print("# Running...")
    args = get_args()
    print("# args:", args)
    pdf_name = main(args)
    print(f'# Done! Made pdf: {pdf_name}')
    # Display timing


#ffmpeg -framerate 20 -pattern_type glob -i '*.png'  -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4