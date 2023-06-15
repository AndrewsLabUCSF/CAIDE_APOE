# Compares b/se using z score method - doi: 10.1038/s41380-019-0596-9
fzscore <- function(b1, b2, se1, se2){(b1 - b2)/sqrt(se1^2 + se2^2)}


farrer_apoe.raw <- read_csv("~/Dropbox/Research/UCSF/ABA/Farrer1997apoex.csv")

farrer_apoe <- farrer_apoe.raw %>%
  mutate(
    b = log(or),
    se = (log(uci) - log(lci)) / 3.92, 
    reference = ifelse(or == 1, TRUE, FALSE),
    apoe = fct_relevel(apoe, "e2/e2", "e2/e3", "e3/e3", "e2/e4", "e3/e4", "e4/e4"), 
    apoe = fct_recode(apoe, "ε2/ε2" = "e2/e2", "ε2/ε3"= "e2/e3", "ε3/ε3\n(reference)" = "e3/e3",
                      "ε2/ε4" = "e2/e4", "ε3/ε4" = "e3/e4", "ε4/ε4" = "e4/e4")
  )


farrer_apoe_diff <- farrer_apoe %>%
  select(race, apoe, b, se) %>%
  pivot_wider(
    id_cols = apoe, names_from = race, values_from = c(b, se)
  ) %>%
  mutate(
     # White vs Black
     diff_nhw_bla = fzscore(b_white, b_black, se_white, se_black), 
     p_diff_nhw_bla = 2*pnorm(q=abs(diff_nhw_bla), lower.tail=FALSE),
     # White vs Hispanic
     diff_nhw_his = fzscore(b_white, b_hispanic, se_white, se_hispanic), 
     p_diff_nhw_his = 2*pnorm(q=abs(diff_nhw_his), lower.tail=FALSE),
     # White vs Japanese
     diff_nhw_jpa = fzscore(b_white, b_japanese, se_white, se_japanese), 
     p_diff_nhw_jpa = 2*pnorm(q=abs(diff_nhw_jpa), lower.tail=FALSE),

     # Black vs Hispanic
     diff_bla_his = fzscore(b_black, b_hispanic, se_black, se_hispanic), 
     p_diff_bla_his = 2*pnorm(q=abs(diff_bla_his), lower.tail=FALSE),
     # Black vs Japanese
     diff_bla_jpa = fzscore(b_black, b_japanese, se_black, se_japanese), 
     p_diff_bla_jpa = 2*pnorm(q=abs(diff_bla_jpa), lower.tail=FALSE),
     
     # Hispanic vs Japanese
     diff_his_jpa = fzscore(b_hispanic, b_japanese, se_hispanic, se_japanese), 
     p_diff_his_jpa = 2*pnorm(q=abs(diff_his_jpa), lower.tail=FALSE)
  ) %>%
  select(apoe, starts_with(c('diff', 'p_diff')))

farrer_apoe_diff_long <- full_join(
  farrer_apoe_diff %>%
    select(apoe, starts_with('diff')) %>%
    pivot_longer(
      cols = starts_with('diff'),
      names_to = 'contrast', 
      values_to = 'z', 
      names_prefix = "diff_",
    ),
  
  farrer_apoe_diff %>%
    select(apoe, starts_with('p_')) %>%
    pivot_longer(
      cols = starts_with('p_'),
      names_to = 'contrast', 
      values_to = 'p', 
      names_prefix = "p_diff_",
    )
) %>%
  separate(contrast, sep = "_", into = c("race1", "race2")) %>%
  mutate(
    stars = case_when(
      p > 0.05 ~ "", 
      between(p,  0.01, 0.05) ~ "*",
      between(p, 0.001, 0.01) ~ "**",
      p < 0.001 ~ "***"
    )
  )

## Data frame for plotting pvalues
farrer_df_pval <- filter(farrer_apoe_diff_long, p < 0.05) %>%
  rename(group1 = race1, group2 = race2) %>%
  mutate(y.position = c(12, 4.5, 4.5, 7.5, 4.5, 12, 4.5, 7.5, 4.5), 
         group1 = case_when(
           group1 == "nhw" ~ "white", 
           group1 == "bla" ~ "black", 
           group1 == "jpa" ~ "japanese", 
           group1 == "his" ~ "hispanic"
         ), 
         group2 = case_when(
           group2 == "nhw" ~ "white", 
           group2 == "bla" ~ "black", 
           group2 == "jpa" ~ "japanese", 
           group2 == "his" ~ "hispanic"
         ),
         race = NA, 
         apoe = fct_relevel(apoe, "e2/e2", "e2/e3", "e3/e3", "e2/e4", "e3/e4", "e4/e4"), 
         apoe = fct_recode(apoe, "ε2/ε2" = "e2/e2", "ε2/ε3"= "e2/e3", "ε3/ε3\n(reference)" = "e3/e3",
                           "ε2/ε4" = "e2/e4", "ε3/ε4" = "e3/e4", "ε4/ε4" = "e4/e4")
  )

theme.size = 14
geom.text.size = (theme.size - 2) * 0.36


farrer_apoe.p <- ggplot(farrer_apoe, aes(y = or, x = race, color = race)) + 
  facet_grid(cols = vars(apoe)) + 
  geom_hline(yintercept = 1, linetype = 2, color = 'grey60') + 
  geom_point(position = position_dodge(width = 0.5)) + 
  scale_y_continuous("Odds Ratio", trans = 'log', 
                     breaks = c(0.25, 0.5, 1, 2, 4, 8, 16, 32), 
                     labels = c(0.25, 0.5, 1, 2, 4, 8, 16, 32)
  ) + 
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")) + 
  geom_errorbar(aes(x = race, ymin = lci, ymax = uci), position = position_dodge(width = 0.5), width = 0) + 
  # geom_vline(xintercept=seq(0.5, 11.5, 1),color="grey90") +
  ggtitle("APOE-related risk for AD") +
  theme_bw() + 
  theme(
    text = element_text(size = theme.size), 
    legend.text = element_text(size = theme.size),
    plot.title = element_text(size=theme.size),
    axis.title.x = element_blank(), 
    axis.text.x = element_text(angle = 315, vjust = 0.5, hjust = 0, size=8),
    axis.text.y = element_text(size=10),
    legend.position = 'bottom', 
    strip.background = element_blank(), 
    panel.grid.minor.y = element_blank(), 
    panel.grid.major.x = element_blank(), 
    legend.title = element_blank(), 
    # panel.border = element_blank(),
  ) + 
  ggprism::add_pvalue(farrer_df_pval, bracket.nudge.y = 3, label = "{stars}", colour = 'black', tip.length = 0.01, bracket.size = 0.2) 


ggsave("results/figures/farrer_apoe_race.png", plot = farrer_apoe.p, width = 5, height = 4, units = 'in')















