# Check the overlapping between the 5% of the nodes with highest
# "Betweenness" CI values and the "Degree"-based hubs

# -------- Samples with Cancer -------------
#Betweenness centrality
betweenness_C <- sort(betweenness(gC), decreasing = T)

#Find the Hubs (top 5%):
hubs_C_betw <- betweenness_C[1:floor(0.05 * length(betweenness_C))]
#and their names:
namesHUBS_C_betw <- names(hubs_C_betw)


#Let's see which "Betweenness"-based hubs are also "Degree"-based hubs
intersect(namesHUBS_C, namesHUBS_C_betw)






# -------- Samples withOUT Cancer -------------
#Betweenness centrality
betweenness_N <- sort(betweenness(gN), decreasing = T)

#Find the Hubs (top 5%):
hubs_N_betw <- betweenness_N[1:floor(0.05 * length(betweenness_N))]
#and their names:
namesHUBS_N_betw <- names(hubs_N_betw)


#Let's see which "Betweenness"-based hubs are also "Degree"-based hubs
intersect(namesHUBS_N, namesHUBS_N_betw)
