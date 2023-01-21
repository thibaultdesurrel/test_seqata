export MipSolver, solve
# unexport : buildMipModel

"""
    MipSolver

Résoud le problème frontale de manière exacte en PLNE.

Exploite l'hypothèse de sujet Seqata, à savoir
que la fonction de coût est linéaire en deux morceaux.

Description du modèle utilisé

### Les données

```
- P (planes) : ensemble des avions p (d'attributs p.lb, p.target, p.ub, ...)
  i est l'indice de p dans les variables de décision
- cost[p,t] : coût de pénalité de l'avion p s'il atterrit à la date t
```

### Les variables de décision

```
x[i in 1:n] : 
date d'atterrissage effective de l'avion i

costs[i in 1:n] : 
coût de chaque avion i compte tenu de sa date d'atterrissage effective

b[i,j] : boolean vrai si l'avion i atterrit avant l'avion j
```
### L'objectif

```
Minimiser total_cost = sum(costs[i]  for i in 1:n
```

### Les contraintes

    AU BOULOT !
    ...

"""
mutable struct MipSolver
    inst::Instance
    bestsol::Solution     # meilleure Solution rencontrée
    # Les attribut spécifiques au modèle
    model::Model  # Le modèle MIP
    x         # vecteur des variables d'atterrissage
    b         # vecteur des variables de précédence (before)
    cost      # variable du cout de la solution
    costs     # variables du cout de chaque avion

    # Le constructeur
    function MipSolver(inst::Instance)
        ln2("MipSolver : constructeur avec $(Args.get("external_lp_solver"))")
        this = new()
        this.inst = inst
        this.bestsol = Solution(inst)
        solver = Args.get("external_lp_solver")
        # Création et configuration du modèle selon le solveur interne sélectionné
        this.model = new_lp_model(mode=:mip, log_level=1)
        return this
    end
end

# Création du modèle frontal pour le problème Seqata (linéaire en deux morceaux)
#
function buildMipModel(sv::MipSolver)

    error("\n\nMéthode buildMipModel(MipSolver) non implantée : AU BOULOT :-)\n\n")
    # ...

end

# Résoud le problème complet : calcul l'ordre, le timing (dates d'atterrissage)
# et le cout total de la solution optimale du problème (coûts linéaires)
# Cette fonction crée le modèle PLNE, le résoud puis met l'objet solution à
# jour.
#
function solve!(sv::MipSolver)
    ln2("BEGIN solve!(MipSolver)")
    ln2( "="^60 )

    lg2("Construction du modèle ($(ms())) ... ")
    buildMipModel(sv)
    ln2("fait ($(ms())).")

    lg2("Lancement de la résolution ($(ms())) ... ")
    optimize!(sv.model)
    ln2("fait ($(ms())).")

    lg2("Test de validité de la solution ($(ms())) ... ")
    if JuMP.termination_status(sv.model) != MOI.OPTIMAL
        print("ERREUR : pas de solution pour :\n    ")
        @show JuMP.termination_status(sv.model)
        println(to_s(sv.bestsol))
        exit(1)
    end
    ln2("fait ($(ms())).")

    lg2("Exploitation des résultats (mise à jour de la solution)($(ms())) ... ")
    # Extraction des valeurs entières des variables d'atterrissage

    # Il reste maintenant à mettre à jour notre objet solution à partir du
    # résultat de la résolution MIP
    for (i, p) in enumerate(sv.inst.planes)
        sv.bestsol.planes[i] = p
        sv.bestsol.x[i] = round(Int, value(sv.x[p.id]))  # car value() retourne Float64 !
        sv.bestsol.costs[i] = value(sv.costs[p.id])
    end
    sv.bestsol.cost = round(value(sv.cost), digits=Args.get(:cost_precision))

    # On trie juste la solution par date d'atterrissage croissante des avions
    # pour améliorer la présentation de la solution
    sort!(sv.bestsol)
    ln2("fait. ($(ms()))")

    ln2("END solve!(MipSolver)")
end
