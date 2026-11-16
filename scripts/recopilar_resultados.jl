# MA4703 — Recolector de resultados A–F (Julia 1.10+).
# Ejecutar DESDE LA RAÍZ DEL REPOSITORIO:
#   julia --project=. scripts/recopilar_resultados.jl
# Para omitir la ejecución de tests:
#   julia --project=. scripts/recopilar_resultados.jl --sin-tests
# No ejecuta simulaciones, optimizaciones, ni renderiza videos.
# No elimina ni modifica resultados existentes.

using Dates
using JLD2
using Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT = joinpath(ROOT, "docs", "paquete_consolidacion")
const DATA = joinpath(ROOT, "data", "processed")
const FIGS = joinpath(ROOT, "docs", "figures")
const VIDEOS = joinpath(ROOT, "docs", "videos")
const RUN_TESTS = !("--sin-tests" in ARGS)

# Que el paquete no se incluya accidentalmente a sí mismo en la copia.
mkpath(OUT)
for name in ("tests", "datos", "curvas", "figuras", "documentacion")
    mkpath(joinpath(OUT, name))
end

rel(path) = replace(relpath(path, ROOT), '\\' => '/')
shortrepr(x; maxlen=180) = begin
    s = replace(repr(x), '\n' => ' ', '\r' => ' ')
    length(s) <= maxlen ? s : first(s, maxlen) * "…"
end

function value_summary(x)
    if x === nothing || x isa Number || x isa Bool || x isa AbstractString || x isa Symbol
        return shortrepr(x)
    elseif x isa NamedTuple
        return shortrepr(x; maxlen=450)
    elseif x isa AbstractArray
        info = "$(typeof(x)); tamaño=$(size(x)); elementos=$(length(x))"
        if length(x) == 0
            return info
        end
        if eltype(x) <: Real && length(x) <= 2_000_000
            lo, hi, finite_count = Inf, -Inf, 0
            for v in x
                if isfinite(v)
                    finite_count += 1
                    lo = min(lo, Float64(v))
                    hi = max(hi, Float64(v))
                end
            end
            info *= "; finitos=$(finite_count)/$(length(x))"
            if finite_count > 0
                info *= "; min=$(lo); max=$(hi)"
            end
        end
        if length(x) <= 12
            info *= "; contenido=" * shortrepr(x; maxlen=250)
        end
        return info
    else
        return "tipo=$(typeof(x)); " * shortrepr(x; maxlen=150)
    end
end

function test_stage(name, test_file)
    logfile = joinpath(OUT, "tests", "test_$(name).txt")
    cmd = `$(Base.julia_cmd()) --project=$(ROOT) $(test_file)`
    ok = false
    elapsed = @elapsed begin
        open(logfile, "w") do io
            try
                run(pipeline(cmd, stdout=io, stderr=io))
                ok = true
            catch err
                println(io, "\nERROR DEL PROCESO: ", sprint(showerror, err))
            end
        end
    end
    println("  Etapa $name: ", ok ? "APROBADA" : "FALLÓ", " (", round(elapsed,digits=1), " s)")
    return (name=name, ok=ok, elapsed=elapsed, filename=logfile)
end

function summarize_jld2(path)
    outfile = joinpath(OUT, "datos", replace(basename(path), r"\.jld2$" => "_resumen.txt"))
    keys_found = String[]
    errors = String[]
    open(outfile, "w") do report
        println(report, "Archivo: ", rel(path))
        println(report, "Tamaño: ", filesize(path), " bytes")
        println(report, "--- Campos y resúmenes (no se re-simula nada) ---")
        try
            JLD2.jldopen(path, "r") do data
                ks = sort!(String.(collect(keys(data))))
                append!(keys_found, ks)
                for k in ks
                    try
                        x = data[k]
                        println(report, k, " = ", value_summary(x))
                    catch err
                        msg = "No se pudo leer $k: " * sprint(showerror,err)
                        push!(errors, msg)
                        println(report, msg)
                    end
                end
            end
        catch err
            msg = "No se pudo abrir el archivo: " * sprint(showerror,err)
            push!(errors, msg)
            println(report, msg)
        end
    end
    return (file=path, summary=outfile, fields=keys_found, errors=errors)
end

# Exportación ligera de curvas escalares para comparar resultados sin Julia/JLD2.
const TIME_GROUPS = [
    ("times", ["error_E", "lyapunov_V", "phase_r", "phase_r2", "lambda_samples"]),
    ("sync_hjb_time", ["sync_hjb_E", "sync_hjb_lambda"]),
    ("sync_constant_time", ["sync_constant_E"]),
    ("desync_hjb_t", ["desync_hjb_r2"]),
]

function export_curves(path)
    produced = String[]
    try
        JLD2.jldopen(path, "r") do d
            for (timekey, cols) in TIME_GROUPS
                timekey in keys(d) || continue
                t = d[timekey]
                t isa AbstractVector{<:Real} || continue
                xs = Pair{String,Any}[]
                for key in cols
                    if key in keys(d)
                        v = d[key]
                        if v isa AbstractVector{<:Real} && length(v) == length(t)
                            push!(xs, key => v)
                        end
                    end
                end
                isempty(xs) && continue
                dest = joinpath(OUT, "curvas", replace(basename(path), ".jld2"=>"") * "_$(timekey).csv")
                stride = max(1, cld(length(t), 501))
                open(dest, "w") do io
                    println(io, join([timekey; first.(xs)], ','))
                    for i in 1:stride:length(t)
                        println(io, join([string(t[i]); [string(pair.second[i]) for pair in xs]], ','))
                    end
                    if length(t)>0 && ((length(t)-1) % stride != 0)
                        i=length(t)
                        println(io, join([string(t[i]); [string(pair.second[i]) for pair in xs]], ','))
                    end
                end
                push!(produced, dest)
            end
        end
    catch err
        println("  Aviso: no se pudieron exportar algunas curvas de $(basename(path)): ", sprint(showerror,err))
    end
    return produced
end

function list_files(directory, ending)
    isdir(directory) || return String[]
    paths = String[]
    for (d, _, files) in walkdir(directory)
        for f in files
            endswith(lowercase(f),lowercase(ending)) && push!(paths, joinpath(d,f))
        end
    end
    return sort!(paths)
end

function human_bytes(bytes)
    bytes >= 1024^3 && return @sprintf("%.1f GB", bytes/1024^3)
    bytes >= 1024^2 && return @sprintf("%.1f MB", bytes/1024^2)
    bytes >= 1024 && return @sprintf("%.1f KB", bytes/1024)
    return "$(bytes) B"
end

function inventory_file(path, files)
    open(path, "w") do io
        println(io,"archivo_relativo,tamano_bytes")
        for f in files
            println(io, "\"", replace(rel(f), '"'=>"\"\""), "\",", filesize(f))
        end
    end
end

function copy_figures(figs; max_total=50*1024^2, max_file=10*1024^2)
    # Se priorizan E y F; si alcanzan, se incluyen A–D. Excluimos archivos demasiado grandes.
    prioritized = sort(figs, by=f->(startswith(basename(f),"F_") ? 0 :
                        startswith(basename(f),"E_") ? 1 : 2, basename(f)))
    selected=String[]
    used=0
    for file in prioritized
        n=filesize(file)
        if n <= max_file && used+n <= max_total
            cp(file,joinpath(OUT,"figuras",basename(file));force=true)
            push!(selected, file)
            used += n
        end
    end
    return selected
end

println("\n=== CONSOLIDACIÓN PROYECTO MA4703 ===")
println("Repositorio: ", ROOT)
println("No se ejecutarán optimizaciones ni simulaciones costosas.")

# Evita que JLD2 tenga problemas al leer structs del proyecto en resultados antiguos.
mainmodule=joinpath(ROOT,"src","FitzHughNagumoControl.jl")
if isfile(mainmodule)
    try
        include(mainmodule)
    catch err
        println("Aviso: no se pudo cargar el módulo para lectura de tipos antiguos: ", sprint(showerror,err))
    end
end

testfiles = [
    ("A", "runtests.jl"),
    ("B", "test_stage_B.jl"),
    ("C", "test_stage_C.jl"),
    ("D", "test_stage_D.jl"),
    ("E", "test_stage_E.jl"),
    ("F", "test_stage_F.jl"),
]
results=NamedTuple[]
if RUN_TESTS
    println("\nEjecutando tests A–F (sin modificar datos ni figuras)...")
    for (name, file) in testfiles
        path = joinpath(ROOT,"tests",file)
        if isfile(path)
            push!(results, test_stage(name,path))
        else
            println("  No encontrado: ",rel(path))
        end
    end
else
    println("\nTests omitidos por --sin-tests.")
end

println("\nLeyendo resultados JLD2...")
files = list_files(DATA,".jld2")
summaries = NamedTuple[]
curves = String[]
for path in files
    push!(summaries,summarize_jld2(path))
    append!(curves,export_curves(path))
    println("  ", basename(path))
end

figs=list_files(FIGS,".png")
vids=vcat(list_files(VIDEOS,".mp4"), list_files(VIDEOS,".gif"))
selected=copy_figures(figs)
inventory_file(joinpath(OUT,"inventario_figuras.csv"),figs)
inventory_file(joinpath(OUT,"inventario_videos.csv"),vids)

# Documentación ya creada por el equipo, sin sobrescribirla.
for file in (joinpath(ROOT,"docs","RESULTADOS_GLOBALES.md"),
             joinpath(ROOT,"docs","ETAPA_E.md"),
             joinpath(ROOT,"docs","ETAPA_F.md"))
    if isfile(file)
        cp(file,joinpath(OUT,"documentacion",basename(file));force=true)
    end
end

reportfile=joinpath(OUT,"RESUMEN_AUTOMATICO.md")
open(reportfile,"w") do io
    println(io,"# MA4703 — Resumen automático de resultados A–F")
    println(io,"\nGenerado: ", Dates.format(now(),"yyyy-mm-dd HH:MM:SS"))
    println(io,"\nEl reporte recoge archivos existentes. No demuestra optimalidad de PMP ni convergencia de HJB.")
    println(io,"\n## Validación")
    if isempty(results)
        println(io,"\nNo se ejecutaron tests en esta recopilación.")
    else
        println(io,"\n| Etapa | Resultado | Tiempo (s) | Log |\n|---|---|---:|---|")
        for r in results
            println(io, "| ", r.name, " | ",
        r.ok ? "Aprobada" : "**Falló**",
        " | ", round(r.elapsed, digits=1),
        " | tests/", basename(r.filename), " |")
        end
    end
    println(io,"\n## Archivos JLD2")
    println(io,"\n| Archivo | Tamaño | Campos | Errores de lectura |\n|---|---:|---:|---:|")
    for s in summaries
        println(io,"| ",rel(s.file)," | ",human_bytes(filesize(s.file))," | ",length(s.fields)," | ",length(s.errors)," |")
    end
    println(io,"\nTodos los campos y rangos de datos están en `datos/`. Las series temporales reducidas están en `curvas/`.")
    println(io,"\n## Gráficos y videos")
    println(io,"\n- Figuras PNG encontradas: **",length(figs),"**; incluidas en el paquete: **",length(selected),"**.")
    println(io,"\n- Videos/GIF encontrados: **",length(vids),"**. No se adjuntan automáticamente para evitar un ZIP gigantesco.")
    println(io,"\n- Archivos de inventario: `inventario_figuras.csv` e `inventario_videos.csv`.")
    println(io,"\n## Precauciones de interpretación")
    println(io,"\n- La condición `E(T) ≤ δ` debe verificarse independientemente en cada candidato de sincronización.")
    println(io,"\n- Un residuo PMP grande impide presentar un candidato como extremal verificado.")
    println(io,"\n- Una HJB fuera de cobertura no suministra una política numérica validada.")
    println(io,"\n- Solo comparar costos cuando `N`, parámetros, condiciones iniciales, restricciones y pesos son iguales.")
end

println("\nResumen generado: ",reportfile)
println("Figuras encontradas: ",length(figs),"; incluidas: ",length(selected))
println("Videos encontrados: ",length(vids)," (inventariados, no copiados)")

zipfile=joinpath(ROOT,"COPTI_Resultados_Para_Revision.zip")
if Sys.iswindows()
    try
        # Generamos ZIP con PowerShell, preinstalado en Windows; el comodín preserva contenido de OUT.
        psquote(s) = "'" * replace(s,"'"=>"''") * "'"
        command = "Compress-Archive -Path " * psquote(joinpath(OUT,"*")) *
                  " -DestinationPath " * psquote(zipfile) * " -Force"
        run(Cmd(["powershell.exe", "-NoProfile", "-NonInteractive", "-Command", command]))
        println("\nZIP CREADO: ",zipfile)
    catch err
        println("\nNo se pudo crear ZIP automáticamente: ",sprint(showerror,err))
        println("Comprímelo manualmente desde la carpeta: ",OUT)
    end
else
    println("\nComprímelo manualmente desde la carpeta: ",OUT)
end
println("\nSube a ChatGPT el ZIP de resultados. Puedes enviar videos MP4 por separado si deseas una revisión visual.")
