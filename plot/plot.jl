# PTRA: Patient Trajectory Analysis Library
# Copyright (c) 2022 imec vzw.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version, and Additional Terms
# (see below).

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public
# License and Additional Terms along with this program. If not, see
# <https://github.com/ExaScience/ptra/blob/master/LICENSE.txt>.

using DataFrames, Statistics, Gadfly, Cairo, Fontconfig, CSV

# Support functions

# Function to compute mean age + stdev, returning a tuple that lists mean age, min age, and max age from a vector 
mean_age(a, m=mean(a), s=std(a)) = (MeanAge=m, MinAge=m-s, MaxAge=m+s)

# Function to plot mean age + std
function plot_cluster_mean(df, title)
    lvls = []
    for i in 0:length(df.CID)-1
        lvls = push!(lvls, i)
    end
    p = plot(df, x=:CID, y=:MeanAge, color=:CID,
    ymin=:MinAge, ymax=:MaxAge, Geom.errorbar, Stat.dodge, Geom.bar(position=:dodge), #Order of Geom on the plot matters for output
    Scale.color_discrete(levels=lvls), #has to be the labels as they occur in the dataframe
    Guide.colorkey(title="clusters"),
    Scale.x_discrete,
    Theme(bar_spacing=0mm, stroke_color=c->"black"),
    Guide.yticks(ticks=[0,10,20,30,40,50,60,70,80,90,100]),
    Guide.Title(title))
    return p
end

# Function to count nr of males/females/ratio
function sex_distr_unique(df::DataFrame)
    fctr = 0
    mctr = 0
    pid_seen = Dict{Int, Bool}();
    for i in 1:length(df.Sex)
        sex = df.Sex[i]
        pid = df.PID[i]
        if !haskey(pid_seen, pid)
            pid_seen[pid] = true
            if sex == "M"
                mctr = mctr + 1
            else
                fctr = fctr + 1
            end
        end
    end
    return (Females=fctr, Males=mctr, Ratio=mctr/fctr)
end

function sex_distr_unique(v1, v2)
    fctr = 0
    mctr = 0
    pid_seen = Dict{Int, Bool}();
    for i in 1:length(v1)
        pid = v1[i]
        sex = v2[i]
        if !haskey(pid_seen, pid)
            pid_seen[pid] = true
            if sex == "M"
                mctr = mctr + 1
            else
                fctr = fctr + 1
            end
        end
    end
    return (Females=fctr, Males=mctr)
end

# Scale male/female counts with regard to total nr of males/females.
function vscale(v, t)
    nv = []
    for i in 1:length(v)
        nv = push!(nv, 100 / t * v[i])
    end
    return nv
end

# Plot sex distribution per cluster
function plot_cluster_sexdistr(df, title)
    plot(df, y=:Sex, x=:Freq, color=:Sex, ygroup=:CID,
        Geom.subplot_grid(Geom.bar(position=:stack, orientation=:horizontal),
            Guide.ylabel(orientation=:vertical)),
        Guide.colorkey(title="Sex"),
        Guide.ylabel("Cluster"), Guide.xlabel("Frequency(%)"),
        Guide.title(title))
end

# Transform sex distribution information into a dataframe for plotting.
function transform_sex_table(df)
    s = []
    f = []
    c = []
    for r in eachrow(df)
        c = push!(c, r.CID)
        s = push!(s, "Female")
        f = push!(f, r."Females(%)")
        c = push!(c, r.CID)
        s = push!(s, "Male")
        f = push!(f, r."Males(%)")
    end
    return DataFrame(CID=c, Sex=s, Freq=f)
end

# Main plot function
function plot_ptra_clusters(pfile, cfile)
    # Parse csv data
    pdf = DataFrame(CSV.File(pfile))
    cdf = DataFrame(CSV.File(cfile))
    println("Patients: ", nrow(pdf))
    println("Clusters: ", nrow(cdf))
    # Join person table and cluster table based on PID
    pcdf = innerjoin(pdf, cdf, on = :PID)
    println("Joined data: ", nrow(pcdf))
    # Filter out entries with unknown age of event of interest
    age_test(age) = age != -1
    pcdf = filter(:AgeEOI => age_test, pcdf)
    println("Filtered Joined data: ", nrow(pcdf))
    pdf = filter(:AgeEOI => age_test, pdf)
    println("Filtered patients: ", nrow(pdf))
    # Group per cluster
    gdf = groupby(pcdf, :CID)

    # Extract per cluster the mean age + std
    agedf = combine(gdf, :Age=>mean_age=>AsTable)
    # Extract per cluster the mean age of the event of interest + std
    ageEOIdf = combine(gdf, :AgeEOI=>mean_age=>AsTable)

    # Plot mean age and mean age of EOI + standard deviations
    set_default_plot_size(21cm, 8cm)
    println("Plotting Mean Age")
    p1 = plot_cluster_mean(agedf, "Mean Age")
    println("Plotting Mean Age EOI")
    p2 = plot_cluster_mean(ageEOIdf, "Mean Age EOI")

    p_sex = sex_distr_unique(pdf)
    println("Male/female distribution: ", p_sex)
    println("Computing make/female distribution per cluster")
    p_unique_sex = combine(gdf, [:PID,:Sex]=>sex_distr_unique=>AsTable) #Counting unique males/females per cluster

    p_unique_sex_pct = select(p_unique_sex, :CID=>:CID, 
                             :Females=>(v->vscale(v, p_sex.Females))=>"Females(%)", 
                             :Males=>(v->vscale(v, p_sex.Males))=>"Males(%)")

    sex_table = transform_sex_table(p_unique_sex_pct)
    println("Plotting sex distribution per cluster")
    p3 = plot_cluster_sexdistr(sex_table, "Sex distribution per cluster (%)")
    return (p1,p2,p3)
end

# For nice plotting of p3: set_default_plot_size(21cm, 25cm)
# For nice plotting of p1, p2: set_default_plot_size(21cm, 8cm)