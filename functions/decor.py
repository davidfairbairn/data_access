# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 11:32:16 2023

@author: 39351
"""
import cartopy as cart
import matplotlib.pyplot as plt 

def decor(ax):
    
    param = ['ocean', 'land']
    for i in param:
        facecolor = i
        if i == 'ocean':
            facecolor='water'
        feature = cart.feature.NaturalEarthFeature(category='physical', name=i, scale='10m',
                                    edgecolor='none', facecolor=cart.feature.COLORS[facecolor], zorder=-1)
        ax.add_feature(feature)
    
    # coastlines
    ax.coastlines(resolution='10m',zorder=-1)
    return ax