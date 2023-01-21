export LpTimingSolver, symbol, solve!

"""
    LpTimingSolver

Résoud du Sous-Problème de Timing par programmation linéaire.

Ce solveur résoud le sous-problème de timing consistant à trouver les dates
optimales d'atterrissage des avions à ordre fixé.
Par rapport aux autres solvers (e.g DescentSolver, AnnealingSolver, ...), il
ne contient pas d'attribut bestsol

VERSION -t lp4 renommé en Lp : modèle diam version sans réoptimisation (coldrestart)
  - modèle diam simplifié par rapport à lp3 : sans réoptimisation (coldrestart)
  - pas de réoptimisation : on recrée le modèle à chaque nouvelle permu
    d'avion dans le solver
  - seules les contraintes de séparation nécessaires à la permu sont créées
  - gestion de l'option --assume_trineq true|false (true par défaut, cf lp1)
  - contraintes de coût simple (diam) : une seule variable de coût par avion
    plus un contrainte par segment :
       cost[i] >= tout_segment[i]
"""
mutable struct LpTimingSolver
    inst::Instance
    # Les attributs spécifiques au modèle
    model::Model  # Le modèle MIP
    x         # vecteur des variables d'atterrissage
    cost      # variable du coût de la solution
    costs     # variables du coût de chaque avion

    nb_calls::Int    # POUR FAIRE VOS MESURES DE PERFORMANCE !
    nb_infeasable::Int

    # Le constructeur
    function LpTimingSolver(inst::Instance)
        this = new()

        this.inst = inst

        # Création et configuration du modèle selon le solveur externe sélectionné
        this.model = new_lp_model() # SERA REGÉNÉRÉ DANS CHAQUE solve!()

        this.nb_calls = 0
        this.nb_infeasable = 0

        return this
    end
end

# Permettre de retrouver le nom de notre XxxxTimingSolver à partir de l'objet
function symbol(sv::LpTimingSolver)
    return :lp
end

function solve!(sv::LpTimingSolver, sol::Solution)

    #error("\n\nMéthode solve!(sv::LpTimingSolver, ...) non implantée: AU BOULOT :-)\n\n")

    sv.nb_calls += 1
    inst = sv.inst
    planes = inst.planes
    n = inst.nb_planes

    # Grandeurs dont on a besoin :
    # Vecteur t de l'ordre des avions donné en argument
    t = [p.id for p in sol.planes]
    # Vecteur E des débuts des fenêtres temporelles d'atterrissage pour tous les avions
    E = [p.lb for p in inst.planes]
    # Vecteur L des fins des fenêtres temporelles d'atterrissage pour tous les avions
    L = [p.ub for p in inst.planes]
    # Matrice S des intervalles de temps nécessaires entre deux atterrissages consécutifs en fonction des types d'avion
    S = inst.sep_mat
    # Vecteur T des heures d'atterrissage souhaitées par avion
    T = [p.target for p in inst.planes]
    # Vecteur K des types d'avion de l'instance
    K = [p.kind for p in inst.planes]
    # Vecteur ep des pénalités unitaires d'avance
    ep = [p.ep for p in inst.planes]
    # Vecteur tp des pénalités unitaires de retard
    tp = [p.tp for p in inst.planes]

    #
    # 1. Création du modèle spécifiquement pour cet ordre d'avion de cette solution
    #
    model = new_lp_model()

    # À COMPLÉTER : variables ? contraintes ? ...
    # ...
    @variable(model, x[1:n], Int)
    # variables définies dans l'énoncé : heures d'atterrissage des avions
    # La variable x sera affectée à l'objet solution, autrement dit, elle doit être suivre le même ordre (et non pas celui de l'instance)
    sv.x = x

    @variable(model, adv[1:n], Int)
    # variable correspondant à l'avance de l'avion : on cherchera à la rendre égale à (T_i - x_i)^+
    # De même, même ordre que la solution

    @variable(model, lat[1:n], Int)
    # variable correspondant au retard de l'avion : on cherchera à la rendre égale à (x_i - T_i)^+
    # Même ordre que la solution

    @expression(model, costs[i in 1:n], ep[i]*adv[i] + tp[i]*lat[i])
    # Cette variable dérivée modélise le coût de pénalité de chaque avion
    sv.costs = costs

    @expression(model, cost, sum(costs))
    # Cette variable modélise l'objectif : le coût de pénalité total de l'horairisation
    sv.cost = cost

    @constraint(model, [i in 1:n], x[i] >= E[i])
    # On s'assure que l'avion t[i] arrive après la borne inférieure de sa TW

    @constraint(model, [i in 1:n], x[i] <= L[i])
    # On s'assure que l'avion t[i] arrive avant la borne inférieure de sa TW

    @constraint(model, [i in 1:n-1], x[t[i+1]] >= x[t[i]] + S[K[t[i]], K[t[i+1]]])
    # On s'assure que les délais de séparation sont respectés ; puisque l'inégalité triangulaire est vérifiée et que 
    # l'ordre des avions est connu, il suffit de les vérifier sur les avions consécutifs.

    @constraint(model, [i in 1:n], adv[i] >= 0)
    @constraint(model, [i in 1:n], adv[i] >= T[i] - x[i])
    # Variable majorant l'avance (positive) de l'avion i par rapport à sa target ; 
    # à l'optimal, elle y est égale (problème de minimisation)

    @constraint(model, [i in 1:n], lat[i] >= 0)
    @constraint(model, [i in 1:n], lat[i] >= x[i] - T[i])
    # Variable majorant le retard (positive) de l'avion t[i] par rapport à sa target ;
    # à l'optimal, elle y est égale (problème de minimisation)

    @objective(model, Min, sum(ep[i]*adv[i] + tp[i]*lat[i] for i in 1:n))
    # Objectif : minimisation du coût total de pénalité des avions

    # 2. résolution du problème à permu d'avion fixée
    #
    JuMP.optimize!(model)

    # 3. Test de la validité du résultat et mise à jour de la solution
    if JuMP.termination_status(model) == MOI.OPTIMAL
        # tout va bien, on peut exploiter le résultat

        # 4. Extraction des valeurs des variables d'atterrissage
        
        sv.x = [value(x[i]) for i in 1:n]
        sv.costs = [ep[i]*value(adv[i]) + tp[i]*value(lat[i]) for i in 1:n]
        sv.cost = sum(sv.costs)

        # ATTENTION : les tableaux x et costs sont dans l'ordre de
        # l'instance et non pas de la solution !
        for (i, p) in enumerate(sol.planes)
            sol.x[i] = round(Int, value(sv.x[p.id]))
        end
        # Mise à jour des coûts (par avion et global) de la solution à partir
        # des dates d'atterrissage. On sous-traite cette mise à jour à la méthode
        # update_costs de la "classe" Solution
        # update_costs!(sol) # diam : serait couteux car recalcule les pénalités
        update_costs!(sol, add_viol_penality=false) # diam => gain de 6% sur alp13

    else
        # La solution du solver est invalide : on utilise le placement au plus
        # tôt de façon à disposer malgré tout d'un coût pénalisé afin de pouvoir
        # continuer la recherche heuristique de solutions.
        sv.nb_infeasable += 1
        solve_to_earliest!(sol)
    end
end
