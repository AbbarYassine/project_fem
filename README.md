					Projet réalisé par Yassine ABBAR dans le cadre du cours de maillage et éléments finis MAIN5 (Polytech Sorbonne)


Comment lancer le projet : (exemple avec disques.geo)
  1) gmsh disques.geo -2 #pour générer le maillage 
  OU placer votre maillage directement dans le dossier du projet si vous l'avez déjà
  2) python eqHelmholz.py meshfile {tag_(interieur} {tag_bord_interieu} {tag_bord_exterieur} #pour lancer
  		ex1 : python eqHelmholtz.py disques.msh 10 11 12 
  		ex2 : python eqHelmholtz.py sous_marin.msh 10 11 12 
  3) paraview visualisation_alpha=3.14_k=0.8.vtu #pour visualiser le résultat sous paraview
  OU ouvrir l'application paraview, charger le fichier puis appuyer sur apply

  Fichiers:
  - Fichier eqHelmholtz.py : résolution du problème Helmholtz 
  - Fichier routinesfem.py: les fonctions python : calcul de matrices fem, assemblage, résolution système ... 
  - Fichiers .msh : données maillages gmsh
  - Fichier ReadMeshData.py : lecture des fichiers .msh
  - Fichier createParaview : création des format .vtu (UnstructuredGrid) avec les données : solution u_sol + maillage 