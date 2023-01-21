export SteepestSolver
export finished, solve!, sample_two_shifts
export record_bestsol, get_stats

"""
SteepestSolver

Résoud le problème global par la statégie de descente profonde.

"""
mutable struct SteepestSolver
    inst::Instance
    nb_test::Int          # Nombre total de voisins testé
    nb_move::Int          # nombre de voisins acceptés

    duration::Float64     # durée réelle (mesurée) de l'exécution
    durationmax::Float64  # durée max de l'exécution (--duration)
    starttime::Float64    # heure de début d'une résolution

    cursol::Solution      # Solution courante
    bestsol::Solution     # meilleure Solution rencontrée
    testsol::Solution     # nouvelle solution potentielle

    bestiter::Int
    do_save_bestsol::Bool
    SteepestSolver() = new() # Constructeur par défaut
end

function SteepestSolver(inst::Instance; startsol::Union{Nothing,Solution} = nothing)
    ln3("Début constructeur de SteepestSolver")

    this = SteepestSolver()
    this.inst = inst
    this.nb_test = 0
    this.nb_move = 0

    this.durationmax = 1.0 * 366 * 24 * 3600   # soit 1 année par défaut !
    this.duration = 0.0 # juste pour initialisation
    this.starttime = 0.0 # juste pour initialisation

    if startsol == nothing
        # Pas de solution initiale => on en crée une
        this.cursol = Solution(inst)
        if Args.get(:presort) == :none
            # initial_sort!(this.cursol, presort=:shuffle) # proto
            initial_sort!(this.cursol, presort = :target) # diam
        else
            initial_sort!(this.cursol, presort = Args.get(:presort))
        end
    else
        this.cursol = startsol
        if lg2()
            println("Dans SteepestSolver : this.cursol = this.opts[:startsol] ")
            println("this.cursol", to_s(this.cursol))
        end
    end

    this.bestsol = Solution(this.cursol)
    this.testsol = Solution(this.cursol)
    this.do_save_bestsol = true
    return this
end

# Retourne true ssi l'état justifie l'arrêt de l'algorithme
#
function finished(sv::SteepestSolver)
    sv.duration = time_ns() / 1_000_000_000 - sv.starttime
    too_long = sv.duration >= sv.durationmax
    other = false
    stop = too_long || other
    if stop
        if lg1()
            println("\nSTOP car :")
            println("     sv.duration=$(sv.duration)")
            println("     sv.durationmax=$(sv.durationmax)")
            println("     (à compléter par vos autres critères d'arrêt)")
            println(get_stats(sv))
        end
        return true
    else
        return false
    end
end

function solve!(
    sv::SteepestSolver;
    startsol::Union{Nothing,Solution} = nothing,
    durationmax::Float64 = 0.0,
    nbh::String = "s1"
)
    ln2("BEGIN solve!(SteepestSolver)")

    println("\nSteepestSolver:solve : ")
    println("   à compléter en s'inspirant de descentsolver")
    println("   AU BOULOT ! ")
    exit(1)

    while !finished(sv)
        sv.nb_test += 1
    end # fin while !finished

    ln2("END solve!(SteepestSolver)")
end

function record_bestsol(sv::SteepestSolver; movemsg = "")
    copy!(sv.bestsol, sv.cursol)
    println("À COMPLÉTER POUR SEQATA !")
end
function get_stats(sv::SteepestSolver)
    txt = """
    ==Etat de l'objet SteepestSolver==
    (à compléter !)
    """
end

# END TYPE SteepestSolver
