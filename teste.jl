# Ativar o ambiente de pacotes específico
import Pkg; Pkg.activate("/home/m3g/Documentos/polyisobutileno/raio_de_giro")
using CSV, DataFrames, Plots, KernelDensity, LaTeXStrings, Glob, Colors

# Definir a paleta de cores
colors = parse.(Colorant, ["#1f77b4ff", "#ff7f0eff", "#2ca02cff", "#d62728ff", "#9467bdff", 
                           "#8c564bff", "#e377c2ff", "#7f7f7fff", "#bcbd22ff", "#17becfff"])

# Função para carregar os dados do arquivo .dat fornecido
function load_data(file_path)
    data = CSV.read(file_path, DataFrame, delim=' ', ignorerepeated=true, header=false)
    if ncol(data) == 1
        lines = readlines(file_path)  # Lê todas as linhas do arquivo
        lines = filter(l -> !startswith(l, "#"), lines)  # Remove cabeçalhos ou linhas sem dados
        data_split = [split(line, r"\s+") for line in lines]  # Divide as linhas por espaços
        data = DataFrame(data_split, :auto)  # Cria o DataFrame com colunas automáticas
        for col in names(data)
            data[!, col] = parse.(Float64, data[!, col])  # Converte as colunas para Float64
        end
    else
        # Apenas converta para Float64 se os dados não forem numéricos
        for col in names(data)
            if eltype(data[!, col]) <: AbstractString  # Se a coluna for do tipo string, converte
                data[!, col] = parse.(Float64, data[!, col])
            end
        end
    end
    return data
end


# Função para extrair a concentração do nome do arquivo
function extract_concentration(file_name)
    parts = split(file_name, '_')  # Separa o nome por "_"
    for part in parts
        if occursin(r"^\d+$", part)  # Procura por partes que sejam números inteiros (concentrações)
            return parse(Int, part)  # Retorna a concentração como inteiro
        end
    end
    return 0  # Retorna 0 se não encontrar a concentração
end

# Função para definir a legenda com base na concentração de coc e oct
function define_legend(concentration)
    if concentration == 0
        return "coc 100%"
    elseif concentration == 20
        return "coc 80%"
    elseif concentration == 40
        return "coc 60%"
    elseif concentration == 60
        return "coc 40%"
    elseif concentration == 80
        return "coc 20%"
    elseif concentration == 100
        return "oct 100%"
    else
        return "unknown"
    end
end

# Listar todos os arquivos .dat no diretório atual
file_paths = Glob.glob("*.dat")

# Dicionário para agrupar os arquivos por parâmetro extraído do nome
file_groups = Dict{String, Vector{String}}()

# Loop para organizar os arquivos no dicionário `file_groups`
for file_path in file_paths
    base_name = splitext(basename(file_path))[1]  # Extrai o nome base do arquivo (sem extensão)
    param = "unknown"  # Define todos como "unknown" inicialmente (parâmetro genérico)
    
    if !haskey(file_groups, param)
        file_groups[param] = []  # Inicializa o grupo se não existir
    end
    push!(file_groups[param], file_path)  # Adiciona o arquivo ao grupo "unknown"
end

# Função para reorganizar o dicionário `file_groups` baseado nas concentrações extraídas
function reorganize_by_concentration(file_groups)
    new_groups = Dict{Int, Vector{String}}()  # Novo dicionário, onde as chaves serão concentrações inteiras
    for (param, files) in file_groups
        for file in files
            concentration = extract_concentration(basename(file))  # Extrai a concentração do nome do arquivo
            if !haskey(new_groups, concentration)
                new_groups[concentration] = []  # Cria um novo vetor vazio, se não existir
            end
            push!(new_groups[concentration], file)  # Adiciona o caminho do arquivo ao grupo correspondente
        end
    end
    return new_groups  # Retorna o novo dicionário reorganizado
end

# Executa a reorganização dos arquivos por concentração
reorganized_groups = reorganize_by_concentration(file_groups)

# Função para plotar gráficos com as configurações fornecidas
    function plot_for_param(files::Vector{String}, colors::Vector{<:Colorant})
        data_list = []  # Lista para armazenar os dados carregados

        # Iterar sobre os arquivos e carregar os dados
        for (i, file_path) in enumerate(files)
            data = load_data(file_path)  # Carrega os dados do arquivo
            base_name = splitext(basename(file_path))[1]  # Extrai o nome base
            concentration = extract_concentration(base_name)  # Extrai a concentração do nome
            legend_label = define_legend(concentration)  # Define a legenda com base na concentração
            push!(data_list, (data, legend_label, concentration, colors[i]))  # Adiciona dado, legenda, concentração e cor à lista
        end

        # Ordenar os dados por concentração (ascendente)
        data_list = sort(data_list, by = x -> x[3])  # Ordenar com base no valor numérico de `concentration`

        # Criação do gráfico combinado
        p1 = plot(legend=:topleft, dpi=500, title="Radius of Gyration vs Frames - Curvas Sobrepostas")
        p2 = plot(legend=:topleft, dpi=500, title="Median Distance from Polymer Center of Mass to Monomer Center")
        plot_font = "Computer Modern"
        default(fontfamily=plot_font, framestyle=:box, label=nothing, grid=false)

        # Iterar sobre os dados e adicionar cada curva ao mesmo gráfico com as configurações fornecidas
        # Iterar sobre os dados e adicionar cada curva ao mesmo gráfico
    for (i, (data, legend_label, _, color)) in enumerate(data_list)
        # Radius of Gyration vs Frames (gráfico p1)
        plot!(p1, data[:, 1] / 100, data[:, 2],
            label=legend_label,
            line=:solid,
            linecolor=color,
            title="Median Polymer Center Mass - Monomer vs Frames",
            xlabel=L"Frames/100",
            ylabel=L"Distance~CM~Polymer~to~Monomer~CM~(Å)",
            xticks=0.0:5.0:maximum(data[:, 1]/100) + 5,
            yticks=10:5.0:maximum(data[:, 2]) + 5,
            xlim=[0.0, maximum(data[:, 1]/100) + 5],
            ylim=[10.0, maximum(data[:, 2]) + 5],
            framestyle=:box,
            grid=true,
            ticks=true,
            linewidth=1,
            legend=:topright)

        # Salvar o gráfico p1 com a nova curva adicionada
        Plots.savefig(p1, "plot_radius_of_gyration_step_$(i).png")

        # Criando gráfico de KDE para a distribuição do Radius of Gyration (gráfico p2)
        kde = KernelDensity.kde(data[!, 2])
        x = range(minimum(data[:, 2]) - 5, stop=maximum(data[:, 2]) + 5, length=length(data[:, 2]))
        plot!(p2, x, z -> pdf(kde, z),
            label="$(legend_label)",
            line=:solid,
            linecolor=color,
            title="",
            xlabel=L"Mediam~distance~cm~poly~to~mon~cm~(Å)",
            ylabel="Probability density",
            xticks=minimum(data[:, 2])-5:10.0:maximum(data[:, 2]) + 5,
            #xtick_labels = round.(xticks, digits=3),
        

            yticks=0.0:0.02:maximum(pdf(kde, x)) + 1,
            xlim=[minimum(data[:, 2])-5, maximum(data[:, 2]) + 1],
            ylim=[0.0, maximum(pdf(kde, x)) + 0.01],
            framestyle=:box,
            ticks=true,
            linewidth=3,
            legend=:topright)

        # Salvar o gráfico KDE p2 com a nova curva adicionada
        Plots.savefig(p2, "plot_kde_step_$(i).png")
    end       




    # Salvar os gráficos individuais e combinados
    Plots.savefig(p1, "multi_plot_combined_$(split(files[1], "/")[end]).png")
    Plots.savefig(p2, "kde_plot_combined_$(split(files[1], "/")[end]).png")
    display(p1)  # Mostrar o gráfico de Radius of Gyration
    display(p2)  # Mostrar o gráfico de KDE
end

# Iterar sobre os grupos reorganizados e gerar gráficos para cada concentração
for (concentration, files) in reorganized_groups
    if length(files) >= 1  # Verificar se há pelo menos 1 arquivo para a concentração
        # Gerar gráficos individuais e o multiplot combinado
        plot_for_param(files, colors[1:length(files)])  # Passar a paleta de cores para o plot
    end
end