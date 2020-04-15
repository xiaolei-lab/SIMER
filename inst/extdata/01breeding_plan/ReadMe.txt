###NOTE###
A breeding plan contains 4 parts:
   1. Genotyping
     (1) This part should start with "# start genotyping"
     (2) Then it follows a table with four columns which names are
         generation, family_index, within_family_index, and sex
     (3) Generation should be no more than "num.gen"
     (4) Numbers should be separated by ","
     (5) ":" can choose serial number, e.g. "1:3,7" represents "1,2,3,7"
     (6) "all" will choose all of the options in a category
     (7) Sex (1 represents sir and 2 represents dam)
     (8) All words should be separated by "Tab"
   2. Phenotyping
     (1) This part should start with "# start phenotyping"
     (2)~(8) are the same as 1. Genotyping 
   3. Fixed effect
     (1) This part should start with "fixed effect"
     (2) Each row represents a trait
     (3) The names of fixed effects should be separated by "Tab"
   4. Random effect
     (1) This part should start with "random effect"
     (2) Each row represents a trait
     (3) The names of random effects should be separated by "Tab"
