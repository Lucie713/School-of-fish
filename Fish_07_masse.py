import streamlit as st
import pandas as pd
import numpy as np
import sys
from PIL import Image, ImageDraw
import math 
import time
from random import randint, gauss
from constants import canvasheigth, canvaswidth, n_grid, sigma, cut_distance, dt

# le 0 est selon +x, et on est dans le sens inverse trigo 



def sumv(u,v):
    (a,b)=u
    (c,d)=v
    return (a+c,b+d)

def difv(u,v):
    (a,b)=u
    (c,d)=v
    return (a-c,b-d)

def norme2(x,y):
    return math.sqrt(x**2+y**2)


def to_cartesian(r, theta):
    """Returns the x and y coordinates of a point defined 
    by its radius r and angle theta (in degrees)"""
    x = r * math.cos(math.radians(theta))
    y = r * math.sin(math.radians(theta))  #bizarre (on a rajouté le -)
    return x, y

class Fish2:

    """classe d'un poisson"""

    def __init__(self,x,y,D,r,R1,R2,R3,dtheta,ecart_type,theta=0,v=200,ax=0,ay=0):
        self.x=x
        self.y=y
        self.D=D
        self.r=r
        self.theta=theta
        self.v=v
        self.vx, self.vy = to_cartesian(self.v, self.theta)
        self.ax=ax
        self.ay=ay
        self.voisinsC1R=[]
        self.voisnsC2AL=[]
        self.voisinsC3AT=[]
        self.R1=R1
        self.R2=R2
        self.R3=R3
        self.dtheta=dtheta
        self.ecart_type=ecart_type
 
       
  
    def draw(self,draw):
        normev=norme2(self.vx,self.vy)
        # on def le poisson avec trois points : le premier est le point de direction, à une distance r du centre 
        draw.polygon(((self.x+self.r*(self.vx/normev),self.y+self.r*(self.vy/normev)),
                     (self.x+(-self.r/2)*(self.vx/normev)+(-self.r/2)*(-self.vy/normev),self.y+(-self.r/2)*(self.vy/normev)+(-self.r/2)*(self.vx/normev)),
                     (self.x+(-self.r/2)*(self.vx/normev)+(self.r/2)*(-self.vy/normev), self.y+(-self.r/2)*(self.vy/normev)+(self.r/2)*(self.vx/normev))),
                       fill=128, outline=128, width=1)
          
        
    def drawcircles(self, draw):
        def draw_main_and_pseudo_circles(R,color,x,y, width):     #x,y : le centre du cercle 
            draw.ellipse([(x-R, y-R),(x+R, y+R)], outline=color, width=width)       # vrai cercle
            #cercles de l'espace torique 
            if x<R:
                draw.ellipse([(x + canvaswidth - R, y - R),(x + canvaswidth + R, y + R)], outline=color, width=width)
            if canvaswidth-x < R:
                draw.ellipse([(x - canvaswidth - R, y - R),(x - canvaswidth + R, y + R)], outline=color, width=width)
            if y <R :
                draw.ellipse([(x  - R, y + canvasheigth - R),(x  + R, y + canvasheigth + R)], outline=color, width=width)
            if canvasheigth - y < R:
                draw.ellipse([(x  - R, y - canvasheigth - R),(x  + R, y - canvasheigth + R)], outline=color, width=width)

        if self.voisinsC1R != []: 
            colorC1='red'  
            widthC1=4
        else : 
            colorC1='indianred' 
            widthC1=3
        if self.voisinsC2AL != []: 
            colorC2='aqua'  
            widthC2=4
        else: 
            colorC2='cornflowerblue'
            widthC2=3
        if self.voisinsC3AT != []: 
            colorC3='limegreen'
            widthC3=4
        else : 
            colorC3='seagreen'
            widthC3=3

        draw_main_and_pseudo_circles(self.R1, colorC1, self.x,self.y, widthC1)
        draw_main_and_pseudo_circles(self.R2, colorC2 ,self.x,self.y, widthC2)
        draw_main_and_pseudo_circles(self.R3, colorC3 ,self.x,self.y, widthC3)

    def champ(self,fish_cell_distance):
        """
        Valeur du "champ de densité" au point x et y.  
        La valeur du champ est une fonction gaussienne de la distance entre le point et le poisson.
        """
        return round(np.exp(-fish_cell_distance/2*(sigma)**2),5) # modifier le sigma

    def champ_gaussien(self):       #que pour un poisson bien situé ( )
        for i in range(n_grid):
            for j in range(n_grid):
                fish_cell_distance = self.distance_torique(((i+0.5)*canvaswidth/n_grid,(j+0.5)*canvasheigth/n_grid))
                if fish_cell_distance < cut_distance:
                    self.D[i][j]=self.champ(fish_cell_distance)
                else:
                    self.D[i][j]=0

    def distance(self,obj):
        if type(obj)==tuple:
            d = math.sqrt((self.y-obj[1])**2+(self.x-obj[0])**2)
        elif type(obj)==Fish2:
            d = math.sqrt((self.y-obj.y)**2+(self.x-obj.x)**2)
        return d
    
    def distance_torique(self,obj):
        if type(obj) == tuple:
            dx = self.x - obj[0]
            dy = self.y - obj[1]
            if abs(dx) > canvaswidth/2:
                dx = canvaswidth - dx 
            if abs(dy) > canvasheigth/2:
                dy = canvasheigth - dy 
            return math.sqrt(dx**2+dy**2)

    def move(self,fishes):   
        self.edge() 
        self.voisinsC1R,self.voisinsC2AL,self.voisinsC3AT=self.voisins_dans_3_cercles(fishes)

#définir le vecteur désiré 
        if  self.voisinsC1R != []:
            d=self.repulsion_desiree()

        elif self.voisinsC2AL !=[] or self.voisinsC3AT !=[]:
            d=sumv(self.alignement_desiree(),self.attraction_desiree())
            (a,b)=d
            d=(a/2,b/2)

        else:
            d=(self.vx,self.vy)





        def vector_angle(vector):
            x, y = vector
            r = np.sqrt(x**2 + y**2)
            angle = np.degrees(np.arctan2(y, x))
            return angle if angle >= 0 else angle + 360
            
        angle_d = vector_angle(d)

#tendre vers le vecteur désiré 
        # retourne le plus petit angle (en valeur absolu), entre deux angles 
        def petit_angle(a,b):
            angle=(a-b)%360
            if angle>180:
                angle=360-angle
            return angle
     
        #si la différence entre angle_d et theta est inférieure à dtheta 
        if petit_angle(angle_d,self.theta) < self.dtheta:
            self.theta = angle_d

        #sinon, trouver vers quel côté doit tourner dtheta (+dtheta ou -dtheta ?) 
        else :
            petit_angle_plus=petit_angle(self.theta+self.dtheta,angle_d)
            petit_angle_moins=petit_angle(self.theta-self.dtheta,angle_d)
            if petit_angle_plus<petit_angle_moins:
                self.theta+=self.dtheta
            else: self.theta-=self.dtheta

#actualiser la vitesse 
        self.x += self.vx*dt
        self.y += self.vy*dt
#wiggle
        wiggle=gauss(0,self.ecart_type)      
        self.theta+=wiggle
       

    def edge(self):
        if self.y<0: self.y=canvasheigth
        if self.y>canvasheigth:self.y=0
        if self.x>canvaswidth:self.x=0
        if self.x<0:self.x=canvaswidth
        self.theta = self.theta % 360
        self.vx, self.vy = to_cartesian(self.v, self.theta)

    def cercle_bis(self,voisin,C):
        liste_cardinaux = [(0,- canvasheigth),(0,canvasheigth),(-canvaswidth,0),(canvaswidth,0)]
        for epsilon in liste_cardinaux :
            distance = ((voisin.x - (self.x + epsilon[0]))**2 + (voisin.y - (self.y + epsilon[1]))**2)
            booleen = distance <= C[0]**2 and distance > C[1]**2
            pseudo_voisin = (voisin.x - epsilon[0], voisin.y - epsilon[1])
            if booleen:
                return (booleen, pseudo_voisin)
        return (False,0)   


    def voisins_dans_3_cercles(self,fishes):
        voisinsC1R=[]
        voisinsC2AL=[]
        voisinsC3AT=[]
        C1=(0,self.R1)
        C2=(self.R1,self.R2)
        C3=(self.R2,self.R3)

        for i in range(len(fishes)):
            D=self.distance(fishes[i])
            if D==0: continue

            if D<=self.R1:
                voisinsC1R.append(fishes[i])
            if self.cercle_bis(fishes[i],C1)[0]:   
                voisinsC1R.append(self.cercle_bis(fishes[i],C1)[1])

            if D>self.R1 and D<=self.R2:
                voisinsC2AL.append(fishes[i])
            if self.cercle_bis(fishes[i],C2)[0]:
                voisinsC2AL.append(fishes[i])

            if D>self.R2 and D<=self.R3:
                voisinsC3AT.append(fishes[i])
            if self.cercle_bis(fishes[i],C3)[0]:   
                voisinsC1R.append(self.cercle_bis(fishes[i],C3)[1])
            
            
        return voisinsC1R,voisinsC2AL,voisinsC3AT
    

    def repulsion_desiree(self):
        dr = (0,0)
        for poisson in self.voisinsC1R:
            if type(poisson) == tuple:
                rij_norm =norme2((poisson[0] - self.x), (poisson[1] - self.y))
                rij = ((poisson[0] - self.x)/rij_norm,(poisson[1] - self.y)/rij_norm)
            else:
                rij_norm =norme2((poisson.x - self.x), (poisson.y - self.y))
                rij = ((poisson.x - self.x)/rij_norm,(poisson.y - self.y)/rij_norm)
            dr = sumv(dr,rij)
        return difv((0,0),dr)

    def alignement_desiree(self):
        nself = norme2(self.vx,self.vy)
        dal = (self.vx/nself,self.vy/nself)   # normaliser 
        for poisson in self.voisinsC2AL:
            vj_norm = norme2(poisson.vx, poisson.vy)
            vj = ((poisson.vx)/vj_norm,(poisson.vy)/vj_norm)
            dal= sumv(dal,vj)
        return dal
    
    def attraction_desiree(self):
        dat = (0,0)
        for poisson in self.voisinsC3AT:
            if type(poisson) == tuple:
                rij_norm =norme2((poisson[0] - self.x), (poisson[1] - self.y))
                rij = ((poisson[0] - self.x)/rij_norm,(poisson[1] - self.y)/rij_norm)
            else :
                rij_norm =norme2((poisson.x - self.x), (poisson.y - self.y))
                rij = ((poisson.x - self.x)/rij_norm,(poisson.y - self.y)/rij_norm)
            dat = sumv(dat,rij)
        return dat