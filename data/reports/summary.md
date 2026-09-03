# Near real-time data on the human neutralizing antibody landscape to influenza virus as of the summer of 2026 to inform vaccine-strain selection _(Summary results)_

## Overview
This report summarizes the results of a study that uses sequencing-based neutralization assays to measure titers to influenza viruses with HAs from human seasonal H3N2 and H1N1 strains representative of those circulating in mid-2026 against human sera collected in early to mid 2026.

For background about sequencing-based neutralization assays, see:

  - [Loes et al (2024), *Journal of Virology*](https://doi.org/10.1128/jvi.00689-24)
  - [Kikawa et al (2026), *eLife*](https://doi.org/10.7554/eLife.106811.4)
  - [Detailed experimental protocol on protocols.io](https://www.protocols.io/view/sequencing-based-neutralization-assay-for-influenz-kqdg3xdmpg25)

The full data and computer code are at <{repo_url}>. This report summarizes the key results; a detailed report of all experimental results at a per-plate and per-serum level is [available here](docs:index.html).

## Download titer data
To download the processed data, go to [this subdirectory](repo:results/final_titer_data). Specifically:

 - [human_titers.csv](repo:results/final_titer_data/human_titers.csv): QC-ed titers for each virus/serum pair.
 - [human_viruses.csv](repo:results/final_titer_data/human_viruses.csv): detailed information about the viruses for which titers were measured.
 - [human_sera.csv](repo:results/final_titer_data/human_sera.csv): detailed information about sera for which titers were measured; each serum is listed once in this file. See also [human_sera_multicohort.csv](repo:results/final_titer_data/human_sera_multicohort.csv) for another file that provides additional sera assignments to more fine-grained cohorts (each serum may be listed multiple times in this file).
 - [human_titers_summarized_by_virus.csv](repo:results/final_titer_data/human_titers_summarized_by_virus.csv): median and geometric mean titers against each virus.
 - [human_summary.txt](repo:results/final_titer_data/human_summary.txt): report on number of sera and viruses analyzed with detailed breakdowns per sera set.

## Viral HAs tested
We chose a set of naturally occurring HAs to represent the diversity of H3N2 and H1N1 influenza in humans as of the summer of 2026, also including recent vaccine strains for both subtypes.
Overall we included 148 HAs (82 H3N2, 66 H1N1) as listed [here](repo:results/final_titer_data/human_viruses.csv).
For the assays, we used viruses with HA ectodomains from these strains and the rest of the genes from the lab-adapted A/WSN/1933 strain.

The HAs we chose largely covered the current diversity of human seasonal influenza, as quantified by the fact (shown in the two figures below) that the HA1 proteins of most H3N2 and H1N1 strains sequenced over the last year are either identical or within one amino-acid mutation of an HA in our library:

![All available H3N2 HA1 protein sequences from the last two years, with the x-axis showing the collection date and the colors indicating how similar the sequence is to a HA in our library.](figure:non-pipeline_analyses/flu_circulating_frequencies/results/plots/H3N2_HA1_counts.svg)

![All available H1N1 HA1 protein sequences from the last two years, with the x-axis showing the collection date and the colors indicating how similar the sequence is to a HA in our library.](figure:non-pipeline_analyses/flu_circulating_frequencies/results/plots/H1N1_HA1_counts.svg)

For interactive Nextstrain phylogenetic trees showing the strains included in the library colored by their subclade identities, see the following links:

 - [H3N2](https://nextstrain.org/community/jbloomlab/flu-seqneut-2026@main/H3N2?c=subclade)
 - [H1N1](https://nextstrain.org/community/jbloomlab/flu-seqneut-2026@main/H1N1?c=subclade)

## Human sera tested
Overall we tested 285 human sera, all collected between April and July 2026 from individuals of a wide range of ages.
See [here](repo:results/final_titer_data/human_sera.csv) for full details on these sera; briefly:

 - 40 residual human sera from Seattle Children's Hospital in the USA (abbreviated SCH)
 - 87 residual human sera from University of Washington Medical Center in the USA (abbreviated UWMC)
 - 38 human sera from a blood donor biobank maintained by Creative Testing Solutions via a collaboration of Vitalant Research Institute and the American Red Cross (abbreviated CTS)
 - 120 human sera from the Victorian Infectious Diseases Reference Laboratory in Australia (abbreviated VIDRL). These sera are pre- and post-vaccination, including:

   + 20 pre- and 20 post-vaccination sera from adults given the Fluzone (egg-based) vaccine
   + 20 pre- and 20 post-vaccination sera from adults given the Flucelvax (cell-based) vaccine
   + 20 pre- and 20 post-vaccination sera from elderly individuals given the Fluad (egg-based) vaccine

## H3N2 results

### Median and per-serum titers

The simplest overview of the results is in the interactive figure below, which shows the median (points) and interquartile range (shaded area) titers for all sera against all strains.
This figure has many interactive options described in its legend; please use them, and for extensive exploration open the plot in a new tab with the link at the end of the legend.

![Median and interquartile range titer of the sera against recent H3N2 strains. The black points indicate the median and the shaded region indicates the interquartile range; mouse over points for details. The tree below the plot indicates the subclade to which each viral strain belongs. The dashed orange line indicates the titer to subclade K, which is the cell-based strain included in the 2026-2027 Northern Hemisphere vaccine. Below the plot you can click on different cohorts to show titers only for sera in those cohorts, and use the sliders to subset only to sera from individuals of specific ages or only sera with a median titer across the strains shown that is within a specific range.](embed:results/titer_plots/human_H3N2_recent_interquartile_range_subclade.html)

The above plot shows that several variants of subclade K have reduced titers, most notably variants with mutations at site 223 (eg, V223I) or 222 (eg, R222K).
Both of these mutations are in antigenic region D, see [Liu et al (2026), *medRxiv*](https://pubmed.ncbi.nlm.nih.gov/41757210/) for discussion of the reason that they may have such a substantial impact in subclade K.

Another view of the same data is in the figure below, which now shows individual lines for each serum rather than the interquartile range.
This figure is much busier, but if you interactively mouse over the lines you can trace the titers for individual sera and visualize the remarkable variability in both overall and strain-specific titers.

![Median and per-serum titers against recent H3N2 strains. The black points indicate the median, and each line indicates the titers of a different serum (mouse over the lines to see the details for a serum and trace its titers to each strain). The tree below the plot indicates the subclade to which each viral strain belongs. The dashed orange line indicates the titer to subclade K, which is the cell-based strain included in the 2026-2027 Northern Hemisphere vaccine. Below the plot you can click on different cohorts to show titers only for sera in those cohorts, and use the sliders to subset only to sera from individuals of specific ages or only sera with a median titer across the strains shown that is within a specific range.](embed:results/titer_plots/human_H3N2_recent_individual_sera_subclade.html)

### Titers projected on phylogenetic tree
A helpful way to examine the titers is on an interactive Nextstrain phylogenetic tree.

The figure below shows a Nextstrain tree (available [as a standalone link here](https://nextstrain.org/community/jbloomlab/flu-seqneut-2026@main/H3N2)) of the HA proteins in the library colored by the median titer against all sera.
The sidebar provides numerous options to change the coloring, including coloring by HA genotype or other properties. You can also click on strains for more details, and use the *Measurements* panel below the tree to break down the titers in various ways.

![Interactive tree of the H3N2 HA proteins, colored by the median titer against all sera. Use the interactive toolbar at left to visualize the data in many other helpful views, and also examine the *Measurements* panel below the tree.](embed:https://nextstrain.org/community/jbloomlab/flu-seqneut-2026@main/H3N2?onlyPanels&d=tree,measurements){height=700}

Note that you can also use the *Scatter* option to the left of the tree to stratify titers against strains with specific mutations, such as [in this view showing the median titers against strains with different amino-acid identities at site 223](https://nextstrain.org/community/jbloomlab/flu-seqneut-2026@main/H3N2?c=gt-HA1_223&d=tree,measurements&l=scatter&m_display=mean&onlyPanels&scatterX=median_titer_All_sera&scatterY=gt).

### Pre- and post-vaccination titers
For the VIDRL cohort, we have pre- and post-vaccination titers for adults who received an egg-based (Fluzone) or cell-based (Flucelvax) vaccine, as well as elderly individuals who received an egg-based (Fluad) vaccine.
The post-vaccination sera are collected a median of 21 days after vaccination (range 18-25 days).

The figure below shows the pre- and post-vaccination titers for each group to all of the strains:

![Pre- and post-vaccination titers to each strain. The colors indicate pre- versus post-vaccine sera, and the points represent the median titers while the shaded region is the interquartile range.](embed:results/pre_post_titer_plots/VIDRL_vaccination_H3N2_recent_interquartile_range_subclade.html)

Here is another figure that shows the **fold change** in titer against each strain after vaccination:

![Fold change in titer post-vaccination relative to pre-vaccination. Points represent the median titers while the shaded region is the interquartile range.](embed:results/pre_post_titer_plots/VIDRL_vaccination_H3N2_recent_interquartile_range_fold_change_subclade.html)

For plots that show the same data but at the level of individual sera rather than interquartile ranges, see [here](docs:results/pre_post_titer_plots/VIDRL_vaccination_H3N2_recent_individual_sera_subclade.html) and [here](docs:results/pre_post_titer_plots/VIDRL_vaccination_H3N2_recent_individual_sera_fold_change_subclade.html).

### Key mutations on HA structure
The titer data above show an effect of mutations at sites 223 and 222 on the titers to subclade K variants.
The interactive figure below shows those sites on an interactive HA structure in purple alongside other mutations in subclade K relative to the 2025-2026 vaccine strain in red.
You can use the dropdown below the figure to also color the HA by antigenic region (you will see that 223 and 222 are in antigenic region D whereas most prior mutations in subclade K were in regions A and B); you can also further interact with the structure using [all of the options in the Mol\* viewer documented here](https://molstar.org/viewer-docs/).

![Key sites shown on the HA protein structure. Use the dropdown below the structure to show additional views, and note also that the structure is interactive (click the wrench icon to open a full toolbar).](embed:results/prot-struct-viz/H3_prot_struct_viz.html)

## H1N1 results

### Median and per-serum titers

The interactive figure below shows the median and interquartile range titers for all sera against all strains.

![Median and interquartile range titer of the sera against recent H1N1 strains. The black points indicate the median and the shaded region indicates the interquartile range; mouse over points for details. The tree below the plot indicates the subclade to which each viral strain belongs. The dashed orange line indicates the titer to subclade D.3.1, which is the cell-based strain included in the 2026-2027 Northern Hemisphere vaccine. Below the plot you can click on different cohorts to show titers only for sera in those cohorts, and use the sliders to subset only to sera from individuals of specific ages or only sera with a median titer across the strains shown that is within a specific range.](embed:results/titer_plots/human_H1N1_recent_interquartile_range_subclade.html)

The above figure shows that most D.3.1.1 strains have reduced titers relative to D.3.1 (the cell-based strain in the 2026-2027 Northern Hemisphere vaccine), with titers especially reduced for strains with mutations at site 155 (eg, G155E) or to a lesser extent 157 (eg, S157L).

Notably, the impact of G155E is greatly accentuated for individuals 15-25 years of age compared to other individuals; you can see this by using the sliders below the plot to subset just on sera from individuals in that age range.

Below is a busier figure showing every individual serum as a line; by mousing over lines you can again see that G155E most affects a subset of sera.

![Median and per-serum titers against recent H1N1 strains. The black points indicate the median, and each line indicates the titers of a different serum (mouse over the lines to see the details for a serum and trace its titers to each strain). The tree below the plot indicates the subclade to which each viral strain belongs. The dashed orange line indicates the titer to subclade D.3.1, which is the cell-based strain included in the 2026-2027 Northern Hemisphere vaccine. Below the plot you can click on different cohorts to show titers only for sera in those cohorts, and use the sliders to subset only to sera from individuals of specific ages or only sera with a median titer across the strains shown that is within a specific range.](embed:results/titer_plots/human_H1N1_recent_individual_sera_subclade.html)

### Sera stratified by impact of G155E
As noted in the section above, the effect of mutations like G155E and S157L is especially profound for sera from a subset of individuals, disproportionately those aged ~15 to 25 years.

To visualize that more clearly, the interactive figure below stratifies sera by their relative titers to two different strains. The initial view separates sera where the titer to strain `D.3.1.1:G155E` is 2-fold or more lower than the titer to `D.3.1.1` (green) versus sera where it is not (orange).
As can be seen, G155E has a very strong effect for a subset of sera, and these sera are enriched for individuals aged ~15 to 25 years.
You can use the interactive options to also stratify sera by relative titers to other pairs of strains.

![Titers stratified by relative titers to two different strains. Colors represent sera where titers to a comparator strain have a fold-change in titer relative to a reference strain that exceeds a threshold; the comparator strain, reference strain, and threshold are all selectable in the interactive options below the plot. The plot then shows the median (points) and interquartile range (shaded region) for each set of sera; the small plot at right shows the age distribution for the two sets; the text above shows the number of sera in each category.](embed:results/stratified_titer_plots/human_stratified_H1N1_recent_interquartile_range_subclade.html)

Note the [comparable plot for H3N2](docs:results/stratified_titer_plots/human_stratified_H3N2_recent_interquartile_range_subclade.html) does not show similarly dramatic stratification for the largest-antigenic-effect mutations.

### Titers projected on phylogenetic tree
The figure below shows a Nextstrain tree (available [as a standalone link here](https://nextstrain.org/community/jbloomlab/flu-seqneut-2026@main/H1N1)) of the HA proteins in the library colored by the median titer against all sera.

![Interactive tree of the H1N1 HA proteins, colored by the median titer against all sera. Use the interactive toolbar at left to visualize the data in many other helpful views, and also examine the *Measurements* panel below the tree.](embed:https://nextstrain.org/community/jbloomlab/flu-seqneut-2026@main/H1N1?onlyPanels&d=tree,measurements){height=700}

Note that you can also use the *Scatter* option to the left of the tree to stratify titers against strains with specific mutations, such as [in this view showing the median titers against strains with different amino-acid identities at site 155](https://nextstrain.org/community/jbloomlab/flu-seqneut-2026@main/H1N1?c=gt-HA1_155&d=tree,measurements&l=scatter&m_display=mean&onlyPanels&scatterX=median_titer_All_sera&scatterY=gt).

### Pre- and post-vaccination titers
The figure below shows the pre- and post-vaccination titers for each group to all of the strains:

![Pre- and post-vaccination titers to each strain. The colors indicate pre- versus post-vaccine sera, and the points represent the median titers while the shaded region is the interquartile range.](embed:results/pre_post_titer_plots/VIDRL_vaccination_H1N1_recent_interquartile_range_subclade.html)

Here is another figure that shows the **fold change** in titer against each strain after vaccination:

![Fold change in titer post-vaccination relative to pre-vaccination. Points represent the median titers while the shaded region is the interquartile range.](embed:results/pre_post_titer_plots/VIDRL_vaccination_H1N1_recent_interquartile_range_fold_change_subclade.html)

For plots that show the same data but at the level of individual sera rather than interquartile ranges, see [here](docs:results/pre_post_titer_plots/VIDRL_vaccination_H1N1_recent_individual_sera_subclade.html) and [here](docs:results/pre_post_titer_plots/VIDRL_vaccination_H1N1_recent_individual_sera_fold_change_subclade.html).

### Key mutations on HA structure
The titer data above show an effect of mutations at site 155, and to a lesser extent 157, on the titers.
The interactive figure below shows those sites on an interactive HA structure in indigo alongside other mutations in subclade D.3.1.1 relative to the 2026-2027 vaccine strain in red.
You can use the dropdown below the figure to also color the HA by antigenic region (you will see that 155 and 157 are in antigenic region Sa).

![Key sites shown on the HA protein structure. Use the dropdown below the structure to show additional views, and note also that the structure is interactive (click the wrench icon to open a full toolbar).](embed:results/prot-struct-viz/H1_prot_struct_viz.html)

## Contributors
This study was led by Caroline Kikawa, Andrew Butler, John Huddleston, and [Jesse Bloom](https://jbloomlab.org/).

Contributors include Heidi Peck and Ian Barr (Doherty Institute, Australia); Janet Englund and Kirsten Lacombe (Seattle Children's Hospital); Alex Greninger (University of Washington); Michael Busch, Marion Lanteri, Mars Stone,
and Bryan Spencer (Vitalant Research Institute and the American Red Cross); Sam Turner and Derek Smith (University of Cambridge); and Scott Hensley (University of Pennsylvania).

