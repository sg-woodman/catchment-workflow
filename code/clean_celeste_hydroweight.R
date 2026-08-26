library(tidyverse)
library(here)


cm_hydro <- read_csv(here("output/CELESTE/CELESTE_hydroweight.csv"))


cm_hydro <- cm_hydro |>
    pivot_longer(cols = 3:1254, names_to = "vari", values_to = "vals") |>
    mutate(
        vari = str_remove(vari, "harvest_regen_y"),
        vari = str_remove(vari, "^ndvi_"),
        vari = str_remove(vari, "_distwtd")
    )

cm_hydro_hr <- cm_hydro |>
    pivot_longer(cols = 3:1254, names_to = "vari", values_to = "vals") |>
    filter(str_detect(vari, "harvest_regen_y")) |>
    mutate(vari = str_remove(vari, "harvest_regen_y")) |>
    separate(vari, into = c("year", "type", "scheme", "drop"), sep = "_") |>
    select(-drop) |>
    rename(prop = vals) |>
    mutate(
        year = as.numeric(year),
        type = str_replace(type, "other", "control")
    )


cm_hydro_hr |>
    group_by(site) |>
    nest() |>
    ungroup() |>
    sample_n(6) |>
    unnest(cols = c(data)) |>
    filter(scheme %in% c("lumped", "haiflo", "haifls")) |>
    ggplot(aes(x = year, y = prop, colour = type)) +
    geom_point() +
    facet_wrap(~ site + scheme, scales = "free_y")


cm_hydro_ndvi <- cm_hydro |>
    pivot_longer(cols = 3:1254, names_to = "vari", values_to = "vals") |>
    filter(str_detect(vari, "ndvi")) |>
    filter(!str_detect(vari, "trend")) |>
    mutate(
        scheme = case_when(
            # HAiFLO/HAiFLS must be checked before iFLO/iFLS: "HAiFLO" contains
            # "iFLO" as a substring, so the shorter pattern would otherwise match
            # first and silently fold HAiFLO/HAiFLS rows into iFLO/iFLS.
            str_detect(vari, "iEucO") ~ "iEucO",
            str_detect(vari, "HAiFLO") ~ "HAiFLO",
            str_detect(vari, "iFLO") ~ "iFLO",
            str_detect(vari, "iEucS") ~ "iEucS",
            str_detect(vari, "HAiFLS") ~ "HAiFLS",
            str_detect(vari, "iFLS") ~ "iFLS",
            TRUE ~ "lumped"
        ),
        year = as.factor(as.numeric(str_extract(vari, "\\d+")) + 1983),
        stat = str_extract(vari, "[^_]+$")
    ) |>
    select(-vari) |>
    pivot_wider(names_from = stat, values_from = vals)

cm_hydro_ndvi |>
    group_by(site) |>
    nest() |>
    ungroup() |>
    sample_n(6) |>
    unnest(cols = c(data)) |>
    filter(scheme %in% c("lumped", "haiflo", "haifls")) |>
    ggplot(aes(x = year, y = median, colour = version)) +
    geom_point() +
    geom_line(group = 1) +
    facet_wrap(~ site + version, scales = "free_y")
