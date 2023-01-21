export DescentSolver
export finished, solve!, sample_two_shifts
export record_bestsol, get_stats
export reset!

"""
    DescentSolver
    DescentSolver(inst::Instance; startsol::Union{Nothing,Solution} = nothing

Résoud le problème global par statégie de descente avec choix aléatoire du
voisin.

Paramètres :
- inst: l'instance à résoudre
- startsol: la solution initiale. Celle-ci peut aussi être passée à chaque appel
  à la méthode `solve!``

Amériorations possibles :
- revoir la gestion des options du solver (utiliser les params par clé-valeur)
- plus besoin de gérer cursol (car ici, bestsol suffit contrainrement à
  ExploreSolver)
"""
mutable struct DescentSolver
    inst::Instance
    nb_test::Int          # Nombre total de voisins testés
    nb_test_max::Int      # Nombre maxi de voisins testés
    nb_move::Int          # nombre de voisins acceptés (améliorant ou non)
    nb_reject::Int        # nombre de voisins refusés
    nb_cons_reject::Int   # Nombre de refus consécutifs
    nb_cons_reject_max::Int # Nombre maxi de refus consécutifs

    duration::Float64     # durée réelle (mesurée) de l'exécution
    durationmax::Float64  # durée max de l'exécution (--duration)
    starttime::Float64    # heure de début d'une résolution

    cursol::Solution      # Solution courante
    bestsol::Solution     # meilleure Solution rencontrée
    testsol::Solution     # nouvelle solution potentielle

    bestiter::Int
    do_save_bestsol::Bool
    DescentSolver() = new() # Constructeur par défaut
end

function DescentSolver(inst::Instance; startsol::Union{Nothing,Solution} = nothing)
    ln3("Début constructeur de DescentSolver")

    this = DescentSolver()
    this.inst = inst
    this.nb_test = 0
    this.nb_test_max = 10_000_000_000  # infini
    this.nb_move = 0
    this.nb_reject = 0
    this.nb_cons_reject = 0
    this.nb_cons_reject_max = 10_000_000_000 # infini

    this.bestiter = 0

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
            println("Dans DescentSolver : this.cursol = this.opts[:startsol] ")
            println("this.cursol", to_s(this.cursol))
        end
    end

    this.bestsol = Solution(this.cursol)
    this.testsol = Solution(this.cursol)
    this.do_save_bestsol = true
    return this
end

# reset! : à utiliser pour plusieurs solves consécutifs
function reset!(sv::DescentSolver)
    sv.cursol.cost = Inf
    sv.bestsol.cost = Inf
    sv.testsol.cost = Inf

    sv.nb_test = 0
    sv.nb_move = 0
    sv.nb_reject = 0
    sv.nb_cons_reject = 0
    sv.bestiter = 0
    sv.duration = 0.0 # juste pour initialisation
    sv.starttime = 0.0 # juste pour initialisation
end

# Retourne true ssi l'état justifie l'arrêt de l'algorithme
#
function finished(sv::DescentSolver)
    sv.duration = time_ns() / 1_000_000_000 - sv.starttime
    too_long = sv.duration >= sv.durationmax
    too_many_reject = (sv.nb_cons_reject >= sv.nb_cons_reject_max)
    too_many_test = (sv.nb_test_max >= sv.nb_test)
    stop = too_long || too_many_reject
    if stop
        if lg1()
            println("\nSTOP car :")
            println("     sv.nb_cons_reject=$(sv.nb_cons_reject)")
            println("     sv.nb_cons_reject_max=$(sv.nb_cons_reject_max)")
            println("     sv.duration=$(sv.duration)")
            println("     sv.durationmax=$(sv.durationmax)")
            println("     sv.nb_test=$(sv.nb_test)")
            println("     sv.nb_test_max=$(sv.nb_test_max)")
            println(get_stats(sv))
        end
        return true
    else
        return false
    end
end

function solve!(
    sv::DescentSolver;
    nb_cons_reject_max::Int = 0,
    startsol::Union{Nothing,Solution} = nothing,
    durationmax::Float64 = 0.0,
    nbh::String = "s1"
)
    ln2("BEGIN solve!(DescentSolver)")
    if durationmax != 0.0
        sv.durationmax = Float64(durationmax)
    end

    if startsol != nothing
        sv.cursol = startsol
        copy!(sv.bestsol, sv.cursol) # on réinitialise bestsol à cursol
        copy!(sv.testsol, sv.cursol)
        if lg2()
            println("Dans DescentSolver : sv.cursol = sv.opts[:startsol] ")
            println("sv.cursol : ", to_s(sv.cursol))
        end
    else
        # on garde la dernière solution sv.cursol
    end

    sv.starttime = time_ns() / 1_000_000_000
    if nb_cons_reject_max != 0
        sv.nb_cons_reject_max = nb_cons_reject_max
    end

    if lg3()
        println("Début de solve : get_stats(sv)=\n", get_stats(sv))
    end

    ln1("\niter <nb_test> =<nb_move>+<nb_reject> <movedesc> => bestcost=...")
    println("Solution de base :", sv.cursol.planes)
    println("Cout de base :", sv.cursol.cost)

    while !finished(sv)
        sv.nb_test += 1
        #prevcost = sv.cursol.cost
        vois = first(nbh)
        d = parse(Int, last(nbh))
        k = rand(1:d)
        i = rand(1:sv.inst.nb_planes-k)

        if vois == 's'
            swap!(sv.cursol, i, i+k)
        elseif vois == 't'
            shift!(sv.cursol, i, i+k)
        elseif vois == 'p'
            indices1 = Array((i:i+k))
            indices2 = shuffle(indices1)
            permu!(sv.cursol, indices1, indices2)
        elseif vois == 'd'
            u = Base.rand()
            if u <= 0.5
                swap!(sv.cursol, i, i+k)
            else
                shift!(sv.cursol, i, i+k)
            end
        else
            ln2("Voisinnage non geré")
            init_solver(sv.cursol, sv.cursol.timing_algo_solver)
        end

        sv.nb_move += 1
        degrad = sv.cursol.cost - sv.bestsol.cost
        ln4("degrad=$(degrad)")
        if degrad < 0
            # Ce voisin est meilleur : on l'accepte
            lg3("+")
            # mise a jour  de la meilleure solution
            copy!(sv.bestsol, sv.cursol)
            println("Nouveau meilleur score : ", sv.cursol.cost)
        else # degrad < 0
            # Ce voisin est plus mauvais : on le rejette
            sv.nb_reject += 1
            sv.nb_cons_reject += 1
        end

    end # fin while !finished
    ln2("END solve!(DescentSolver)")
end

# Returne un quadruplet d'indices destiné à affectuer deux shifts relativement
# proches
# - ecartmaxin est l'écart maxi au sein d'une paire d'indices
# - ecartmaxout est l'écart maxi entre deux paires d'indices (cumulables)
# - ecartmaxin et ecartmaxout sont imposés dans les bornes de l'instance
# - abs(i2-i1) et abs(i4-i3) sont limités par ecartmaxin
# - i3 est distant au maximum de ecartmaxout du couple (i1,i2)
function sample_two_shifts(sol::Solution; ecartmaxin::Int = 10, ecartmaxout::Int = -1)

    # **PROTO**
    # Version stupide car voisinage trop large !
    i1 = rand(1:sol.inst.nb_planes)
    i2 = rand(1:sol.inst.nb_planes)
    i3 = rand(1:sol.inst.nb_planes)
    i4 = rand(1:sol.inst.nb_planes)
    return (i1, i2, i3, i4)

end

function record_bestsol(sv::DescentSolver; movemsg = "")
    copy!(sv.bestsol, sv.cursol)
    sv.bestiter = sv.nb_test
    if sv.do_save_bestsol
        write(sv.bestsol)
    end
    if lg3()
        print("\niter $(rpad(sv.nb_test, 4))=$(sv.nb_move)+$(sv.nb_reject) ")
        print("$movemsg ")
        print("bestsol=$(to_s(sv.bestsol))")
    elseif lg1()
        print("\niter $(rpad(sv.nb_test, 4))=$(sv.nb_move)+$(sv.nb_reject) ")
        print("$movemsg => bestcost=", sv.cursol.cost)
    end
end
function get_stats(sv::DescentSolver)
    # txt = <<-EOT.gsub /^ {4}/,''
    txt = """
    ==Etat de l'objet DescentSolver==
    sv.nb_test=$(sv.nb_test)
    sv.nb_test_max=$(sv.nb_test_max)
    sv.nb_move=$(sv.nb_move)
    sv.nb_cons_reject=$(sv.nb_cons_reject)
    sv.nb_cons_reject_max=$(sv.nb_cons_reject_max)

    sv.duration=$(sv.duration)
    sv.durationmax=$(sv.durationmax)

    sv.testsol.cost=$(sv.testsol.cost)
    sv.cursol.cost=$(sv.cursol.cost)
    sv.bestsol.cost=$(sv.bestsol.cost)
    sv.bestiter=$(sv.bestiter)
    sv.testsol.solver.nb_infeasable=$(sv.testsol.solver.nb_infeasable)
    """
    txt = replace(txt, r"^ {4}" => "")
end

# END TYPE DescentSolver
