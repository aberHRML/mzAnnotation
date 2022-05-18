#' Calculate the rings plus double bonds equivalent
#' @description Calculate the rings plus double bonds equivalent (RDBE) for molecular formulas.
#' @param element_frequencies table of element frequencies for a set of molecular formulas as returned by `elementFrequencies()`.
#' @param valences named list of element valences
#' @return A vector of RDBE values.
#' @examples 
#' element_frequencies <- elementFrequencies(c('C12H22O11','C12H22NO11'))
#' rdbe(element_frequencies)
#' @importFrom tibble as_tibble
#' @importFrom dplyr distinct
#' @importFrom tidyr replace_na
#' @export

rdbe <- function(element_frequencies,
                 valences = list(C = 4,
                                 H = 1,
                                 N = 3,
                                 O = 2,
                                 P = 3,
                                 S = 4)){
  valences <- valences %>% 
    as_tibble() %>% 
    gather(element,valence) %>% 
    mutate(valence = valence - 2)
  
  element_frequencies <- element_frequencies %>% 
    gather(element,
           frequency,
           -MF) %>% 
    mutate(frequency = replace_na(frequency,0)) %>% 
    left_join(valences,by = 'element')
  
  rdbe_values <- element_frequencies %>% 
    group_by(MF,element) %>% 
    mutate(value = frequency * valence) %>%
    ungroup() %>% 
    group_by(MF) %>% 
    summarise(rdbe = sum(value) %>% 
                {. + 2} %>% 
                {./2}
    )
  element_frequencies %>% 
    select(MF) %>% 
    dplyr::distinct() %>% 
    left_join(rdbe_values,by = 'MF') %>% 
    select(rdbe) %>% 
    deframe()
}

#' LEWIS and SENIOR checks
#' @rdname lewis_senior
#' @description LEWIS molecular formula valence test and SENIOR test for the existence of molecular graphs. 
#' @param element_frequencies table of element frequencies for a set of molecular formulas as returned by `elementFrequencies()`.
#' @param valences named list of element valences
#' @return Boolean vector of check results for each molecular formula.
#' @examples 
#' element_frequencies <- elementFrequencies(c('C12H22O11','C12H22NO11'))
#' lewis(element_frequencies)
#' senior(element_frequencies)
#' @export

lewis <- function(element_frequencies,
                  valences = list(C = 4,
                                  H = 1,
                                  N = 3,
                                  O = 2,
                                  P = 3,
                                  S = 4)){
  tibble(RDBE = rdbe(element_frequencies,
                     valences = valences)) %>% 
    rowwise() %>% 
    mutate(remainder = RDBE %% 1,
           LEWIS =  ifelse(
             all(
               RDBE >= 0,
               remainder != 0.5
             ),
             TRUE,
             FALSE)) %>% 
    .$LEWIS
}

#' @rdname lewis_senior
#' @export

senior <- function(element_frequencies,
                   valences = list(C = 4,
                                   H = 1,
                                   N = 3,
                                   O = 2,
                                   P = 3,
                                   S = 4)){
  valences <- valences %>% 
    as_tibble() %>% 
    gather(element,valence)
  
  mfs <- element_frequencies$MF
  
  element_frequencies <- element_frequencies %>%  
    gather(element,frequency,-MF) %>% 
    mutate(frequency = replace_na(frequency,0)) %>% 
    left_join(valences,by = 'element') %>% 
    mutate(total_valence = frequency * valence) %>% 
    group_by(MF)
  
  element_frequencies %>% 
    summarise(sum_valence = sum(total_valence)) %>% 
    select(MF,sum_valence) %>% 
    left_join(element_frequencies %>%
                filter((valence %% 2) != 0) %>% 
                summarise(odd_valence_total = sum(frequency)),
              by = 'MF') %>% 
    left_join(element_frequencies %>% 
                filter(frequency > 0) %>% 
                summarise(twice_maximum_valence = max(valence) * 2),
              by = 'MF') %>% 
    left_join(element_frequencies %>% 
                summarise(twice_atoms_minus_1 = sum(frequency) %>% 
                            {. * 2 - 1}),
              by = 'MF') %>% 
    rowwise() %>% 
    mutate(SENIOR = ifelse(
      all(
        (sum_valence %% 2) == 0 | (odd_valence_total %% 2) == 0,
        sum_valence >= twice_maximum_valence,
        sum_valence >= twice_atoms_minus_1
      ),
      TRUE,
      FALSE
    )) %>% 
    mutate(MF = factor(MF,levels = mfs)) %>% 
    arrange(MF) %>% 
    .$SENIOR
}

#' Element ratio checks
#' @description Element ratio checks based on rules #4 and #5 of the Seven Golden Rules by Kind et al 2007.
#' @param element_ratios a tibble containing molecular formula elemental ratios
#' @param range the ratio threshold ranges as defined by rules #5 and #5 of the Seven Golden Rules by Kind et al 2007
#' @return A tibble containing results of the element ratio checks.
#' @examples 
#' elementFrequencies(c('H2O','C12H22O11')) %>% 
#'   elementRatios() %>% 
#'   elementRatioCheck()
#' @importFrom purrr flatten_chr
#' @importFrom dplyr group_split
#' @importFrom rlang parse_expr eval_tidy
#' @export

elementRatioCheck <- function(element_ratios, 
                                 range = c('common',
                                           'extended',
                                           'extreme')){
  range_types <- c('common',
                   'extended',
                   'extreme')
  
  range_type <- match.arg(range,
                          choices = range_types)
  
  thresholds <- tibble(range = range_types %>% 
                         map(rep,times = 6) %>% 
                         flatten_chr(),
                       ratio = rep(c(rep('H/C',2),'N/C','O/C','P/C','S/C'),3),
                       operator = c(rep(c('>',rep('<',5)),2),'<',rep('>',5)),
                       threshold = c(0.2,3.1,1.3,1.2,0.3,0.8,
                                     0.1,6,4,3,2,3,
                                     0.1,6,1.3,1.2,0.3,0.8)) %>% 
    filter(range == range_type) %>%
    mutate(name = paste0(ratio,' ',
                         operator,
                         ' ',threshold),
           expr = paste0('element_ratios[["',ratio,'"]] ',
                         operator,
                         ' ',threshold)
    )
  
  ratio_checks <- thresholds %>% 
    rowwise() %>% 
    group_split() %>% 
    map_dfc(~{
      result <- .x$expr %>% 
        parse_expr() %>% 
        eval_tidy()
      
      if (length(result > 0)) result <- result else result <- NA
      
      checks <- tibble(
        !!.x$name := result
      )
      
      return(checks)
    }) %>% 
    bind_cols(select(element_ratios,MF)) %>% 
    select(MF,everything())
  
  return(ratio_checks)
}

#' Element count checks
#' @description Element count checks based on rule #6 of the Seven Golden Rules by Kind et al 2007.
#' @param element_frequencies a tibble containing element frequencies as returned by `elementFrequencies()`
#' @return A tibble containing results of the element count checks.
#' @examples 
#' elementFrequencies(c('H2O','C12H22O11')) %>% 
#'   elementCountCheck()
#' @importFrom rlang .data
#' @export

elementCountCheck <- function(element_frequencies){
  
  check_names <- c(rep('NOPS all >= 1',4),
                   rep('NOP all >= 3',3),
                   rep('OPS all >= 1',3),
                   rep('PSN all >= 1',3),
                   rep('NOS all >= 6',3))
  
  element_checks <- tibble(
    name = check_names,
    element = c('N','O','P','S',
                'N','O','P',
                'O','P','S',
                'P','S','N',
                'N','O','S'),
    operator = '>=',
    count = c(rep(1,4),
              rep(3,3),
              rep(1,3),
              rep(1,3),
              rep(6,3))
  ) %>% 
    dplyr::mutate(check = paste0(.data$element,' ',
                                 .data$operator,' ',
                                 .data$count),
                  expr = paste0('element_frequencies[["',.data$element,'"]]',
                                ' ',.data$operator, ' ',
                                .data$count))
  
  probability_checks <- check_names %>% 
    unique() %>% 
    map_dfc(~{
      element_freq_checks <- element_checks %>% 
        filter(name == .x) %>% 
        rowwise() %>% 
        group_split() %>% 
        map_dfc(~{
          result <- .x$expr %>% 
            parse_expr() %>% 
            eval_tidy()
          
          if (length(result) > 0) result <- result 
          else result <- NA
          
          checks <- tibble(
            !!.x$check := result)
          
        }) 
      
      if (nrow(element_freq_checks) > 0){
        element_freq_checks %>% 
          rowid_to_column(var = 'row') %>% 
          gather(check,result,-row) %>% 
          group_by(row) %>% 
          summarise(!!.x := all(result)) %>% 
          select(-row)
      } else {
        NULL
      }
    }) %>% 
    rowid_to_column(var = 'row') %>% 
    gather(name,probability_check,-row)
  
  heuristics <- element_checks %>% 
    mutate(
      operator = '<',
      count = c(10,20,4,3,
                11,22,6,
                14,3,3,
                3,3,4,
                19,14,8
      ),
      check = paste0(.data$element,' ',
                     .data$operator,' ',
                     .data$count),
      expr = paste0('element_frequencies[["',.data$element,'"]]',
                    ' ',.data$operator, ' ',
                    .data$count)
    )
  
  heuristic_checks <- check_names %>% 
    unique() %>% 
    map_dfc(~{
      heuristic_checks <- heuristics %>% 
        filter(name == .x) %>% 
        rowwise() %>% 
        group_split() %>% 
        map_dfc(~{
          checks <- tibble(
            !!.x$check := .x$expr %>% 
              parse_expr() %>% 
              eval_tidy()
          )    
          
          if (nrow(checks) == 0) NULL
          else checks
        }) 
      
      if (nrow(heuristic_checks) > 0){
        heuristic_checks  %>% 
          rowid_to_column(var = 'row') %>% 
          gather(check,result,-row) %>% 
          group_by(row) %>% 
          summarise(!!.x := all(result)) %>% 
          select(-row)
      } else {
        tibble(!!.x := NA)
      }
    }) %>% 
    rowid_to_column(var = 'row') %>% 
    gather(name,heuristic_check,-row)
  
  check_results <- left_join(probability_checks,
                             heuristic_checks, 
                             by = c("row", "name")) %>% 
    mutate(heuristic_check = replace(.data$heuristic_check,
                                     .data$probability_check == FALSE,
                                     NA) %>% 
             replace(is.na(.data$probability_check),
                     NA)) %>% 
    group_by(row,.data$name) %>% 
    summarise(result = all(.data$heuristic_check),.groups = 'drop') %>% 
    spread(name,result) %>% 
    select(-row) %>% 
    bind_cols(select(element_frequencies,.data$MF)) %>% 
    select(.data$MF,everything())
  
  heuristic_names <- heuristics %>% 
    select(name,check) %>% 
    group_by(name) %>% 
    summarise(label = paste(check,collapse = ', ')) %>% 
    mutate(label = paste(name,label,sep = '; ')) %>% 
    .$label
  
  check_results <- check_results %>% 
    set_colnames(c('MF',heuristic_names))
  
  return(check_results)
}

#' The proportion of C, H and O atoms in molecular formulas
#' @description Calculate the proportion of C, H and O in molecular formulas.
#' @param element_frequencies a tibble containing element frequencies as returned by `elementFrequencies()`
#' @return A tibble contianing the 
#' @examples
#' elementFrequencies(c('H2O','C12H22O11')) %>% 
#'   CHOproportion()
#' @export

CHOproportion <- function(element_frequencies){
  element_totals <- element_frequencies %>% 
    gather(element,
           count,
           -MF) %>% 
    group_by(.data$MF) %>% 
    summarise(total_atoms = sum(.data$count,na.rm = TRUE))
  
  CHO_counts <- element_frequencies %>% 
    gather(element,
           count,
           -MF) %>% 
    filter(.data$element %in% c('C','H','O'))
  
  
  if (nrow(CHO_counts) == 0) {
    CHO_counts <- tibble(
      MF = element_frequencies$MF,
      CHO = 0
    )
    } else {
      CHO_counts <- CHO_counts %>% 
        group_by(MF) %>% 
        summarise(CHO = sum(.data$count,na.rm = TRUE))
    }
   CHO_counts %>% 
    left_join(element_totals,by = 'MF') %>% 
    mutate(`CHO proportion` = CHO/.data$total_atoms)
}

#' Golden rule tests for molecular formlas
#' @description Heuristic tests for moleucla formulas based on the golden rules 2, 4, 5 and 6 from Kind et al 2007.
#' @param MF a vector of molecular formulas
#' @return A tibble containing golden rule heuristic check results.
#' @examples 
#' goldenRules(c('H2O','C12H22O11'))
#' @export

goldenRules <- function(MF){
  element_frequencies <- elementFrequencies(MF)
  element_ratios <- elementRatios(element_frequencies)
  
  golden_rules <- tibble(
    MF = MF,
    LEWIS = lewis(element_frequencies),
    SENIOR = senior(element_frequencies),
  ) %>% 
    left_join(elementRatioCheck(element_ratios),
              by = 'MF') %>% 
    left_join(elementCountCheck(element_frequencies),
              by = 'MF') %>% 
    left_join(CHOproportion(element_frequencies) %>% 
                select(MF,`CHO proportion`),
              by = 'MF')
  
  return(golden_rules)
}

#' Molecular formula plausibility scores
#' @description Percentage plausibility scores based on rules 2, 4, 5 and 6 Kind et al 2007.
#' @param golden_rules a tibble containing golden rule heuristic checks results as from `goldenRules()`
#' @return A tibble containing golden rules plausibility scores.
#' @examples 
#' goldenRules(c('H2O','C12H22O11')) %>% 
#'   goldenRulesScore()
#' @export

goldenRulesScore <- function(golden_rules){
  
  rule_types <- tibble(
    check = colnames(golden_rules)[-1]
  ) %>% 
    mutate(rule = .data$check,
           rule = replace(.data$rule,.data$rule == 'LEWIS','LEWIS and SENIOR'),
           rule = replace(.data$rule,.data$rule == 'SENIOR','LEWIS and SENIOR'),
           rule = replace(.data$rule,grepl('/',.data$rule),'Element ratios'),
           rule = replace(.data$rule,grepl('all',.data$rule),'Element counts'))
  
  golden_rules %>% 
    gather(check,result,-MF) %>% 
    left_join(rule_types, 
              by = "check") %>% 
    mutate(result = replace(result,
                            result == TRUE,
                            1) %>% 
             replace(is.na(result),
                     1)) %>%
    group_by(MF,rule) %>% 
    summarise(score = sum(result)/dplyr::n(),
              .groups = 'drop') %>% 
    spread(rule,score) %>% 
    mutate(`Plausibility (%)` = (`LEWIS and SENIOR` +
                                   `Element ratios` +
                                   `Element counts` +
                                   `CHO proportion`) / 4 * 100) %>% 
    select(MF,
           `LEWIS and SENIOR`,
           `Element ratios`,
           `Element counts`,
           `CHO proportion`,
           `Plausibility (%)`)
}
