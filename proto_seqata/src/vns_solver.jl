export VnsSolver

"""
VnsSolver

Résoud le problème global par méthode à voisinage variable.

"""
mutable struct VnsSolver
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
    VnsSolver() = new() # Constructeur par défaut
end

function VnsSolver(inst::Instance; startsol::Union{Nothing,Solution} = nothing)
    ln3("Début constructeur de VnsSolver")

    this = VnsSolver()
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
            println("Dans VnsSolver : this.cursol = this.opts[:startsol] ")
            println("this.cursol", to_s(this.cursol))
        end
    end

    this.bestsol = Solution(this.cursol)
    this.testsol = Solution(this.cursol)
    this.do_save_bestsol = true
    return this
end


function finished(sv::VnsSolver)
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
    sv::VnsSolver;
    nb_cons_reject_max::Int = 0,
    startsol::Union{Nothing,Solution} = nothing,
    durationmax::Float64 = 0.0,
    nbh::String = "s2"
)
    ln2("BEGIN solve!(VnsSolver)")
    if durationmax != 0.0
        sv.durationmax = Float64(durationmax)
    end

    if startsol != nothing
        sv.cursol = startsol
        copy!(sv.bestsol, sv.cursol) # on réinitialise bestsol à cursol
        copy!(sv.testsol, sv.cursol)
        if lg2()
            println("Dans VnsSolver : sv.cursol = sv.opts[:startsol] ")
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
        n = sv.inst.nb_planes
        k = 1

        # A cette étape, sv.cursol et sv.bestsol contiennent x
        # Le k va correspondre à la taille des voisinage qu'on considère
        # La variable nbh contiendra le type de voisinage ainsi que k_max
        vois = first(nbh)
        k_max = parse(Int, nbh[2:length(nbh)])
        while k <= k_max
            # Etape de perturbation
            # On génére un voisin aléatoire de x dans V_k(x)
            indice_random = Random.rand(1:n-k)
            if vois == 's'
                swap!(sv.cursol, indice_random, indice_random + k)
            elseif vois == 't'
                shift!(sv.cursol, indice_random, indice_random + k)
            else
                error("Voisinnage "*nbh*" non geré")
                #init_solver(sv.cursol, sv.cursol.timing_algo_solver)
            end
            # A cette étape, sv.cursol contient un voisin x' choisit aléatoirement dans V_k(x)
            # Etape de recherche locale
            # On regarde tous les swap/shift de taille k de x' pour voir si on a une solution améliorante
            for i in 1:n-k
                copy!(sv.testsol, sv.cursol)

                if vois == 's'
                    swap!(sv.testsol, i, i + k)
                elseif vois == 't'
                    shift!(sv.testsol, i, i + k)
                else
                    error("Voisinnage "*nbh*" non geré")
                    #init_solver(sv.cursol, sv.cursol.timing_algo_solver)
                end

                if sv.testsol.cost < sv.cursol.cost
                    copy!(sv.cursol, sv.testsol)
                end
            end
            # A cette étape sv.cursol contient l'optimum local ainsi obtenu

            # Etape de deplacement
            if sv.cursol.cost < sv.bestsol.cost
                copy!(sv.bestsol, sv.cursol)
                k = 1
                println("Nouveau meilleur score : ", sv.cursol.cost)
            else
                k += 1
                sv.nb_cons_reject += 1
            end
        end

    end # fin while !finished
    ln2("END solve!(VnsSolver)")
end
# END TYPE VnsSolver
function get_stats(sv::VnsSolver)
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
