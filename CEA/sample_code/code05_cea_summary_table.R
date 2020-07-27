# Make EIC-style table

# pull out ceaf

cea_table = cbind(icer_table[,c("Strategy", "Cost", "Effect", "ICER_conc")], 
                  ceac_table[,2:5])

for(i in c(2, 5:8)){cea_table[,i] = round(cea_table[,i], 2)}
cea_table[,"Cost"] = format(round(cea_table[,"Cost"], 0), nsmall=0, big.mark=",")
cea_table[,"Effect"] = format(round(cea_table[,"Effect"], 0), nsmall=0, big.mark=",")

cea_table = tibble(cea_table) %>% 
  mutate(`0` = as.character(`0`)) %>% 
  mutate(`250` = as.character(`250`)) %>% 
  mutate(`500` = as.character(`500`)) %>% 
  mutate(`1000` = as.character(`1000`)) 


for (i in 1:length(lbl$int_names)){
  tmp1=as.character(cea_table$Strategy[lns[i]:(lns[i+1]-1)])==as.character(ceaf_list[[i]]$int_short[ceaf_list[[i]]$thresholds==0])
  tmp2=as.character(cea_table$Strategy[lns[i]:(lns[i+1]-1)])==as.character(ceaf_list[[i]]$int_short[ceaf_list[[i]]$thresholds==250])
  tmp3=as.character(cea_table$Strategy[lns[i]:(lns[i+1]-1)])==as.character(ceaf_list[[i]]$int_short[ceaf_list[[i]]$thresholds==500])
  tmp4=as.character(cea_table$Strategy[lns[i]:(lns[i+1]-1)])==as.character(ceaf_list[[i]]$int_short[ceaf_list[[i]]$thresholds==1000])
  
  
  cea_table[lns[i]:(lns[i+1]-1),] = cea_table[lns[i]:(lns[i+1]-1),] %>% 
    mutate(`0` = cell_spec(`0`, "latex", background = ifelse(tmp1, "pink", "white"))) %>% 
    mutate(`250` = cell_spec(`250`, "latex", background = ifelse(tmp2, "pink", "white"))) %>% 
    mutate(`500` = cell_spec(`500`, "latex", background = ifelse(tmp3, "pink", "white"))) %>% 
    mutate(`1000` = cell_spec(`1000`, "latex", background = ifelse(tmp4, "pink", "white"))) 
}

cea_table = as.matrix(cea_table[,2:9])
rownames(cea_table) = c(labelstratnames2_yb, rep(labelstratnames2_cod, times=4))

colnames(cea_table) = c("Mean cost difference", "Mean DALYs averted", "ICER",
                        "Minimum cost", "250 USD per DALY averted", "500 USD per DALY averted", "1,000 USD per DALY averted", "Pr EOT by 2030")

# wrap colnames with str_wrap

tb_kb = kable(cea_table, format="latex", 
              escape=T, row.names=T, 
              caption = paste0("Summary of cost-effectiveness, assuming a time horizon of ", hrzn[1], "-", hrzn[length(hrzn)], ". 
                  Cost differences and DALYs averted are relative to the comparator, which is the first strategy listed for each location. 
                  DALYs averted and cost differences are discounted at 3 percent per year in accordance with guidelines. 
                  In the uncertainty analysis (columns 5-8), the probability that a strategy is optimal is shown 
                  (as a proportion of all simulations, accounting for parameter uncertainty).
                  Strategies highlighted in pink are optimal strategies: the strategies for which the mean net monetary 
                  benefit (NMB) is highest. ICER: Incremental Cost-effectiveness Ratio, 
                  DALY: disability adjusted life-years. For an extended discussion of these terms, 
                  see supplement section \\ref{subsec:supp_cea}."),
              label = paste0("ceasummaries-",hrzn[1], "-", hrzn[length(hrzn)]),
              align = c(rep("R{1.1cm}",2), "R{1.4cm}", rep("R{1.25cm}", 4), "R{1.15cm}"))

tb_kb = kableExtra::pack_rows(tb_kb, placenames[1], 1, 2)
for(i in 2:length(foldernames)){tb_kb = kableExtra::pack_rows(tb_kb, placenames[i], (i-1)*4-1, (i-1)*4+2)}

tb_kb = add_header_above(tb_kb, c(" " = 1, "Cost-effectiveness analysis\nwithout uncertainty" = 3, 
                                  "Net benefit (uncertainty) analysis:
                                    \nPr. that a strategy is optimal,
                                    \n(conditional on willingness-to-pay)" = 4, " " = 1),
                         color="white", background = "black", bold=T)

# if I want to add left or right line to header rows
# https://community.rstudio.com/t/add-left-border-to-additional-header-row-in-kableextra/26330

tb_kb = gsub("[t]", "", tb_kb, fixed=T)
tb_kb = gsub("\\begin{table}", "\\begin{table}[h!]", tb_kb, fixed=T)

tb_kb = row_spec(tb_kb, 0, color="white", background = "black", bold=T)
# NOTE: this line might not work as intended. Add white edges? Or underline.

tb_kb = gsub("\\centering", "\\centering \\footnotesize", tb_kb, fixed=T)
tb_kb = gsub("\\end{tabular}", "\\end{tabular} \\normalsize", tb_kb, fixed=T)

tb_kb = gsub("\\multicolumn{6}", "\\multicolumn{9}", tb_kb, fixed=T)
tb_kb = gsub("\\multicolumn{5}", "\\multicolumn{9}", tb_kb, fixed=T)
tb_kb = gsub("\\multicolumn{9}", "\\rowcolor{black!20} \\multicolumn{9}", tb_kb, fixed=T) ## add colors

tb_kb = gsub("\\rowcolor{black!20} \\multicolumn{9}{l}{\\textbf{Yasa Bonga}}","\\arrayrulecolor{black} \\rowcolor{black!20} \\multicolumn{9}{l}{\\textbf{Yasa Bonga}}", tb_kb, fixed=T) ## add colors

tb_kb = gsub("\\textbackslash{}", "\\", tb_kb, fixed=T)
tb_kb = gsub("\\{", "{", tb_kb, fixed=T)
tb_kb = gsub("\\}", "}", tb_kb, fixed=T)

arraycolorsheader = "\\hhline{|>{\\arrayrulecolor{black}}->{\\arrayrulecolor{white}}|---|----|>{\\arrayrulecolor{black}}-}
\\arrayrulecolor{white}"
tb_kb = gsub("\\cline{2-4} \\cline{5-8}", arraycolorsheader, tb_kb, fixed=T)

tb_kb = gsub("\\makecell[c]{Cost", "\\parbox{4.0cm}{\\centering Cost", tb_kb, fixed=T) ## add colors
tb_kb = gsub("\\makecell[c]{Net", "\\parbox{4.4cm}{\\centering Net", tb_kb, fixed=T) ## add colors

write(tb_kb, file=paste0("/Users/Marina/GD_AYA/Sleeping sickness/DRC5 Analysis/cea_sum", hrzn[1], "-", hrzn[length(hrzn)], ".tex"))
write(tb_kb, file=paste0("/Users/Marina/Dropbox/Apps/Overleaf/DRC 5 SZ CEA manuscript/tables/cea_sum", hrzn[1], "-", hrzn[length(hrzn)], ".tex"))
