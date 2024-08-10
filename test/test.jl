using DataFrames
using Plots

# Exemplo de dados fictícios
# Cada linha representa um país, e cada coluna representa a posição nas Olimpíadas daquele ano
data = DataFrame(
    País = ["Estados Unidos", "China", "Rússia", "Reino Unido", "Alemanha", "França", "Itália", "Austrália", "Coreia do Sul", "Japão",
            "Canadá", "Holanda", "Brasil", "Espanha", "Cuba", "Hungria", "Nova Zelândia", "Polônia", "Noruega", "Suécia"],
    `2004` = [1, 2, 3, 4, 5, 6, 8, 10, 9, 11, 12, 13, 16, 14, 17, 19, 18, 20, 15, 7],
    `2008` = [1, 2, 3, 4, 6, 5, 7, 8, 9, 10, 12, 13, 15, 14, 17, 19, 18, 20, 11, 16],
    `2012` = [1, 2, 4, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 19, 17],
    `2016` = [1, 2, 4, 3, 5, 6, 8, 7, 10, 9, 11, 12, 13, 15, 14, 16, 17, 19, 20, 18],
    `2020` = [1, 2, 3, 4, 5, 6, 8, 9, 7, 10, 12, 13, 11, 15, 14, 16, 17, 18, 20, 19],
    `2024` = [1, 3, 2, 4, 5, 6, 7, 8, 9, 10, 12, 13, 11, 14, 15, 16, 18, 19, 20, 17]
)

# Transpor os dados para o formato adequado para plotagem
df_melted = stack(data, Not(:País))

# Renomear colunas para melhor entendimento
rename!(df_melted, Dict(:variable => :Ano, :value => :Posição))

# Plotar os dados
plot(
    df_melted.Ano, df_melted.Posição,
    group = df_melted.País,
    label = df_melted.País,
    xlabel = "Ano",
    ylabel = "Posição",
    title = "Posição dos 20 Melhores Países ao Longo das Últimas 6 Olimpíadas",
    lw = 2,
    legend = :outertopright,
    xrotation = 45,
    xticks = unique(df_melted.Ano),
    yticks = 1:20,
    reverse = :y  # Inverte a ordem para que 1 esteja no topo
)
