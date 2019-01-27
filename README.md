					Projet réalisé dans le cadre du cours de maillage et éléments finis MAIN5 (Polytech Sorbonne)


Comment lancer le projet : (exemple avec disques.geo)
  1) gmsh Maillage/disques.geo -2 #pour le maillage dans le dossier Maillage 
  OU placer votre maillage directement dans le dossier Maillage si vous l'avez déjà
  2) python eqHelmholz.py #pour lancer le fichier sousmarin_simple.msh 
  3) paraview Output_alpha=alphavalue.vtu #pour visualiser le résultat sous paraview
  OU ouvrir l'application paraview, charger le fichier puis appuyer sur apply

  Fichiers:
  - Fichier eqHelmholtz.py : résolution du problème Helmholtz 
  - Fichier routinesfem: les fonctions python : calcul de matrices fem, assemblage, résolution système ... 
  - Fichiers.msh : données maillages gmsh 