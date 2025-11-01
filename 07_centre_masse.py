import streamlit as st
import pandas as pd
import numpy as np
import altair as alt
import sys
from PIL import Image, ImageDraw
import math 
import time
from random import randint, gauss
from Fish_07_masse import Fish2 
from constants import canvasheigth, canvaswidth, n_grid, dt, temps_caracteristique, temps_boucle


def norme2(x,y):
    return math.sqrt(x**2+y**2)

def sumv2(u,v):
    (a,b)=u
    (c,d)=v
    return (a+c,b+d)

def sumv3(u,v):
    (a,b,c)=u
    (d,e,f)=v
    return (a+d,b+e,c+f)

def difv2(u,v):
    (a,b)=u
    (c,d)=v
    return (a-c,b-d)

def produit_vectoriel(u,v):
    (a,b,c)=u
    (d,e,f)=v
    return (b*f-e*c, c*d-f*a, a*e-d*b)

def polarisation(fishes):
    pgroupex=0
    pgroupey=0
    for fish in fishes :
        v=norme2(fish.vx,fish.vy)
        pgroupex+=fish.vx/v
        pgroupey+=fish.vy/v
    return norme2(pgroupex,pgroupey)/Nbpoissons

def liste_poissons_translates(fishes,centre_densite):
    """ Cette fonction renvoie une liste des positions (sous forme de tuple) des poissons
      translatés de façon à ce que le banc soit situé au centre de l'écran. On opère :
       - une translation de tous les poissons de façon à ce que le centre de densité soit
         confondu avec le centre du cadre. 
       - pour les poissons qui sont sortis du cadre : on fait des modulos selon la taille du cadre"""

    poissons_translates=[ 0 for _ in range (Nbpoissons)]

    for i in range(Nbpoissons):
        nouveau_x = fishes[i].x - centre_densite[0] + canvaswidth/2 
        nouveau_y = fish.y - centre_densite[1] + canvasheigth/2

        if nouveau_x > canvaswidth : 
            nouveau_x -= canvaswidth
        if nouveau_x < 0 :
            nouveau_x += canvaswidth
        
        if nouveau_y > canvasheigth :
            nouveau_y -= canvasheigth
        if nouveau_y < 0 :
            nouveau_y += canvasheigth

        poissons_translates[i] = (nouveau_x, nouveau_y)
    return poissons_translates 



def coordonnees_centre_masse(fishes,centre_densite):
    """ Renvoie un couple qui contient le centre de masse de tous les poissons translatés, 
    et le véritable centre de masse."""
    poissons_translates = liste_poissons_translates(fishes,centre_densite)
    centre_intermediaire = (0,0)

    #on calcule classiquement le centre de masse translaté en sommant toutes les positions des poissons translatés.
    for fish in poissons_translates :
        centre_intermediaire = (centre_intermediaire[0] + fish[0], centre_intermediaire[1] + fish[1])


    # normalisation
    centre_intermediaire = (centre_intermediaire[0]/Nbpoissons, centre_intermediaire[1]/Nbpoissons)

    # On récupère le centre de masse réel 
    centre_masse = (centre_intermediaire[0] + centre_densite[0] - canvaswidth/2, 
                    centre_intermediaire[1] + centre_densite[1] - canvasheigth/2)

    return (centre_intermediaire,centre_masse)


def ric(fish, centre_intermediaire):
    """ Renvoie un vecteur normalisé qui correspond à centre_intermediaire - poissons_translate"""
    (a,b) = difv2((fish[0],fish[1]), centre_intermediaire)
    norme_ric = norme2(a,b)
    return (a/norme_ric,b/norme_ric,0)


def moment_angulaire(fishes,centre_densite):
    moment_groupe=(0,0,0)
    centre_intermediaire = coordonnees_centre_masse(fishes,centre_densite)[0]
    poissons_translates = liste_poissons_translates(fishes,centre_densite)
    for i in range (Nbpoissons) :
        v=norme2(fishes[i].vx,fishes[i].vy)
        moment_groupe=sumv3(moment_groupe, produit_vectoriel(ric(poissons_translates[i],centre_intermediaire), (fishes[i].vx/v,fishes[i].vy/v,0)))
    return abs(moment_groupe[2])/Nbpoissons


#affichage
with st.sidebar:
    with st.form(key='Paramètres'):
     Nbp = st.text_input(label='nombre de poissons',value="2")
     vmi = st.text_input(label='vitesse initiale minimum',value="300")
     vma = st.text_input(label='vitesse initiale maximum',value="300")
     dthe = st.text_input(label='dtheta',value="10")
     R11 = st.text_input(label='R1',value="60")
     R22= st.text_input(label='R2',value="600")
     R33 = st.text_input(label='R3',value="900")
     Ecart_type=st.text_input(label='wiggle',value='5')
     les_cercles = st.checkbox('Afficher les cercles',value=False)
     calculer_la_densite = st.checkbox('Calculer la densité',value=False)
     afficher_la_densite = st.checkbox('Afficher la densité',value=False)
     afficher_centre_de_densite = st.checkbox('Afficher le centre de densité',value=False)
     afficher_centre_de_masse = st.checkbox('Afficher le centre de masse',value=False)
     calculer_moment_angulaire_polarisation = st.checkbox('Calculer le moment angulaire et la polarisation',value=False)
     afficher_polarisation_et_moment_angulaire = st.checkbox('Afficher le moment angulaire et la polarisation à la fin de l animation ',value=False)
     afficher_animation = st.checkbox('Afficher l animation', value=False)
     submit_button = st.form_submit_button(label='submit')
   

Nbpoissons=int(Nbp)
vmin=int(vmi)
vmax=int(vma)
dtheta=int(dthe)
R1=int(R11)
R2=int(R22)
R3=int(R33)
ecart_type=int(Ecart_type)
 



# création de la liste de poissons
fishes=[]
while len(fishes)<Nbpoissons:
  
    p=Fish2(x=randint(canvaswidth//2-100,canvaswidth//2+100),y=randint(canvasheigth//2-100,canvasheigth//2+100),D=np.zeros((101,101)),
         r=randint(30,30),theta=randint(0,360),v=randint(vmin,vmax),R1=R1,R2=R2,R3=R3,dtheta=dtheta,ecart_type=ecart_type)
  #  p=Fish2(x=randint(0,canvaswidth),y=randint(0,canvasheigth),D=np.zeros((101,101)),
  #         r=randint(30,30),theta=randint(0,360),v=randint(vmin,vmax),R1=R1,R2=R2,R3=R3,dtheta=dtheta,ecart_type=ecart_type)
    fishes.append(p)

###########
#val = [(1501, 434, 30), (1090, 983, 69), (565, 519, 34), (780, 843, 49), (234, 1961, 277), (1970, 1549, 281), (1594, 1490, 196), (695, 1813, 143), (1300, 1865, 167), (1898, 931, 122), (1136, 1308, 91), (848, 1009, 139), (519, 1026, 251), (1949, 1043, 334), (1850, 1132, 262), (976, 1871, 358), (1917, 706, 313), (14, 364, 327), (413, 408, 79), (1625, 72, 211), (948, 1806, 147), (608, 325, 156), (1589, 1072, 71), (1857, 688, 146), (706, 1697, 72), (1728, 723, 344), (231, 554, 152), (1307, 180, 83), (1198, 1555, 225), (1679, 1443, 260), (1473, 1411, 287), (1267, 1428, 3), (1518, 216, 68), (1111, 107, 293), (710, 1685, 224), (564, 1267, 69), (128, 82, 279), (1667, 430, 7), (18, 1795, 269), (161, 966, 176), (1371, 225, 333), (1342, 423, 159), (1502, 1860, 78), (1503, 1999, 250), (468, 242, 354), (158, 1982, 158), (963, 1155, 64), (555, 605, 297), (1028, 1791, 221), (18, 1154, 230)]

# création de la liste de poissons
#fishes=[]
#i = 0
#while len(fishes)<Nbpoissons:
#  
#    p=Fish2(x = val[i][0], y = val[i][1],D=np.zeros((101,101)),
#           r=randint(30,30),theta = val[i][2],v=randint(vmin,vmax),R1=R1,R2=R2,R3=R3,dtheta=dtheta,ecart_type=ecart_type)
 #   fishes.append(p)
#    i += 1
############

def matrice_de_densite(fishes):  # somme des gaussiennes de tous les poissons
    Dtot = np.zeros((n_grid, n_grid))
    for fish in fishes:
        fish.champ_gaussien()
        Dtot = sumv_taille(Dtot, fish.D)
    return Dtot

def sumv_taille(A, B):  # fonction de somme matricielle pour des matrices de la bonne taille
    n_rows, n_cols = A.shape
    M = np.zeros((n_rows, n_cols))
    for i in range(n_rows):
        for j in range(n_cols):
            M[i][j] = A[i][j] + B[i][j]
    return M

def draw_density(Dt, draw):
    "dessine la densité du groupe"
    n = len(Dt)
    cell_width = canvaswidth // n
    cell_height = canvasheigth // n
    max_density = np.max(Dt)
    
    for i in range(n):
        for j in range(n):
            density = Dt[i][j] / max_density
            intensity = int(255 * density)
            color = (intensity, intensity, intensity)
            draw.rectangle([i * cell_width, j * cell_height, (i + 1) * cell_width, (j + 1) * cell_height], fill=color)

#création du canvas 
canvas=st.empty()

times=[]
polarization=[]
angular_momentum=[]
for i in range(1,temps_boucle):

    for fish in fishes:
        fish.move(fishes)
    
    if calculer_la_densite and i > temps_caracteristique:
        Dt = matrice_de_densite(fishes)
        argmax = np.argmax(Dt)      # que faire si il y a deux max égaux ? 
        centre_densite_ligne =  argmax//n_grid
        centre_densite_colonne = argmax%n_grid
        centre_densite = (centre_densite_ligne*canvaswidth/n_grid,centre_densite_colonne*canvasheigth/n_grid)
        if calculer_moment_angulaire_polarisation:
                angular_momentum.append(moment_angulaire(fishes,centre_densite))
                polarization.append(polarisation(fishes))


    if afficher_animation:
        img = Image.new("RGB", (canvaswidth, canvasheigth), (230, 230, 230))
        draw = ImageDraw.Draw(img)
        if calculer_la_densite and i > temps_caracteristique:
            if afficher_la_densite:
                draw_density(Dt, draw)
            if afficher_centre_de_densite:
                draw.ellipse([(centre_densite[0] - 5, centre_densite[1] -5 ),(centre_densite[0] +5 , centre_densite[1]+5)], outline='black', width=5)
            if afficher_centre_de_masse:
                (a,b) = coordonnees_centre_masse(fishes,centre_densite)[1]
                draw.ellipse([(a-5,b-5),(a+5,b+5)],outline='red', width=5)
        for fish in fishes:
            if les_cercles :
                fish.drawcircles(draw)
            fish.draw(draw)    

       # st.text(f"theta={fish.theta},vx={fish.vx},vy={fish.vy},angle={math.atan(fish.vy/fish.vx)*360/(2*math.pi)}")
        
    #times : le temps à partir du temps de coupure (il faut rajouter le temps caractéristique pour avoir le temps total d'animation)
    if i > temps_caracteristique:    
        times.append(i*dt)
    
    if afficher_animation:
        canvas.image(img)
  
    time.sleep(dt)



# afficher la grille totale Dtotale à la fin du programme 
#df = pd.DataFrame(matrice_de_densite(fishes))
#st.dataframe(df)


# afficher les graphiques de polarisation et de moment angulaire

if afficher_polarisation_et_moment_angulaire:
    if calculer_la_densite and calculer_moment_angulaire_polarisation and temps_boucle > temps_caracteristique:

        if len(times) == len(polarization) == len(angular_momentum):
            data={"time": times, "polarisation":polarization,"moment_angulaire":angular_momentum}
            df = pd.DataFrame(data)  
        
            base = alt.Chart(df).encode(
                alt.X('time').title("Temps (s)")
            )

            line = base.mark_line(stroke='#5276A7', interpolate='monotone').encode(
                alt.Y('polarisation',scale=alt.Scale(domain=[0, 1])).title('Polarisation', titleColor='#5276A7')
            )

            line2 = base.mark_line(stroke='pink', interpolate='monotone').encode(
                alt.Y('moment_angulaire',scale=alt.Scale(domain=[0, 1])).title('Moment angulaire', titleColor='pink')
            )

            combined = alt.layer(line, line2).resolve_scale(
                y='independent'
            )

            st.altair_chart(combined, use_container_width=True)
        else:
            st.text("Les listes 'times', 'polarization' et 'angular_momentum' doivent avoir la même longueur.")
    else :  
        st.text(f"vous avez pas calculé la densite/moment ang/polarisation  : impossible d'afficher les courbes ")

def average(l):
    if not l:
        return 0
    return sum(l) / len(l)
if calculer_la_densite and calculer_moment_angulaire_polarisation and temps_boucle > temps_caracteristique:
    avg_polarization = average(polarization)
    avg_angular_momentum = average(angular_momentum)
    st.text(f"Polarisation moyenne: {avg_polarization}, Moment angulaire moyen: {avg_angular_momentum}")
else :  
    st.text(f"vous avez pas calculé la densite/moment ang/polarisation : impossible de renvoyer la moyenne")

# st.line_chart(data=df)

