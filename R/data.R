#' @title Spambase
#' @docType data
#' @usage data(spambase)
#' @description This is a well known dataset with a binary target obtainable from the UCI machine learning dataset archive. Each row is an e-mail, which is considered to be either spam or not spam. The dataset contains 48 attributes that measure the percentage of times a particular word appears in the email, 6 attributes that measure the percentage of times a particular character appeared in the email, plus three attributes measuring run-lengths of capital letters.
#' @format A data frame with 4,601 rows and 58 variables (1 categorical, 57 continuous). 
#' \describe{
#'   \item{\code{is.spam}}{Is the email considered to be spam? (0=no,1=yes)}
#'   \item{\code{word.freq.make}}{Percentage of times the word 'make' appeared in the e-mail}
#'   \item{\code{word.freq.address}}{Percentage of times the word 'address' appeared in the e-mail}
#'   \item{\code{word.freq.all}}{Percentage of times the word 'all' appeared in the e-mail}
#'   \item{\code{word.freq.3d}}{Percentage of times the word '3d' appeared in the e-mail}
#'   \item{\code{word.freq.our}}{Percentage of times the word 'our' appeared in the e-mail}
#'   \item{\code{word.freq.over}}{Percentage of times the word 'over' appeared in the e-mail}
#'   \item{\code{word.freq.remove}}{Percentage of times the word 'remove' appeared in the e-mail}
#'   \item{\code{word.freq.internet}}{Percentage of times the word 'internet' appeared in the e-mail}
#'   \item{\code{word.freq.order}}{Percentage of times the word 'order' appeared in the e-mail}
#'   \item{\code{word.freq.mail}}{Percentage of times the word 'mail' appeared in the e-mail}
#'   \item{\code{word.freq.receive}}{Percentage of times the word 'receive' appeared in the e-mail}
#'   \item{\code{word.freq.will}}{Percentage of times the word 'will' appeared in the e-mail}
#'   \item{\code{word.freq.people}}{Percentage of times the word 'people' appeared in the e-mail}
#'   \item{\code{word.freq.report}}{Percentage of times the word 'report' appeared in the e-mail}
#'   \item{\code{word.freq.addresses}}{Percentage of times the word 'addresses' appeared in the e-mail}
#'   \item{\code{word.freq.free}}{Percentage of times the word 'free' appeared in the e-mail}
#'   \item{\code{word.freq.business}}{Percentage of times the word 'business' appeared in the e-mail}
#'   \item{\code{word.freq.email}}{Percentage of times the word 'email' appeared in the e-mail}
#'   \item{\code{word.freq.you}}{Percentage of times the word 'you' appeared in the e-mail}
#'   \item{\code{word.freq.credit}}{Percentage of times the word 'credit' appeared in the e-mail}
#'   \item{\code{word.freq.your}}{Percentage of times the word 'your' appeared in the e-mail}
#'   \item{\code{word.freq.font}}{Percentage of times the word 'font' appeared in the e-mail}
#'   \item{\code{word.freq.000}}{Percentage of times the word '000' appeared in the e-mail}
#'   \item{\code{word.freq.money}}{Percentage of times the word 'money' appeared in the e-mail}
#'   \item{\code{word.freq.hp}}{Percentage of times the word 'hp' appeared in the e-mail}
#'   \item{\code{word.freq.hpl}}{Percentage of times the word 'hpl' appeared in the e-mail}
#'   \item{\code{word.freq.george}}{Percentage of times the word 'george' appeared in the e-mail}
#'   \item{\code{word.freq.650}}{Percentage of times the word '650' appeared in the e-mail}
#'   \item{\code{word.freq.lab}}{Percentage of times the word 'lab' appeared in the e-mail}
#'   \item{\code{word.freq.labs}}{Percentage of times the word 'labs' appeared in the e-mail}
#'   \item{\code{word.freq.telnet}}{Percentage of times the word 'telnet' appeared in the e-mail}
#'   \item{\code{word.freq.857}}{Percentage of times the word '857' appeared in the e-mail}
#'   \item{\code{word.freq.data}}{Percentage of times the word 'data' appeared in the e-mail}
#'   \item{\code{word.freq.415}}{Percentage of times the word '415' appeared in the e-mail}
#'   \item{\code{word.freq.85}}{Percentage of times the word '85' appeared in the e-mail}
#'   \item{\code{word.freq.technology}}{Percentage of times the word 'technology' appeared in the e-mail}
#'   \item{\code{word.freq.1999}}{Percentage of times the word '1999' appeared in the e-mail}
#'   \item{\code{word.freq.parts}}{Percentage of times the word 'parts' appeared in the e-mail}
#'   \item{\code{word.freq.pm}}{Percentage of times the word 'pm' appeared in the e-mail}
#'   \item{\code{word.freq.direct}}{Percentage of times the word 'direct' appeared in the e-mail}
#'   \item{\code{word.freq.cs}}{Percentage of times the word 'cs' appeared in the e-mail}
#'   \item{\code{word.freq.meeting}}{Percentage of times the word 'meeting' appeared in the e-mail}
#'   \item{\code{word.freq.original}}{Percentage of times the word 'original' appeared in the e-mail}
#'   \item{\code{word.freq.project}}{Percentage of times the word 'project' appeared in the e-mail}
#'   \item{\code{word.freq.re}}{Percentage of times the word 're' appeared in the e-mail}
#'   \item{\code{word.freq.edu}}{Percentage of times the word 'edu' appeared in the e-mail}
#'   \item{\code{word.freq.table}}{Percentage of times the word 'table' appeared in the e-mail}
#'   \item{\code{word.freq.conference}}{Percentage of times the word 'conference' appeared in the e-mail}
#'   \item{\code{char.freq.;}}{Percentage of times the character ';' appeared in the e-mail}
#'   \item{\code{char.freq.(}}{Percentage of times the character '(' appeared in the e-mail}
#'   \item{\code{char.freq.[}}{Percentage of times the character '[' appeared in the e-mail}
#'   \item{\code{char.freq.!}}{Percentage of times the character '!' appeared in the e-mail}
#'   \item{\code{char.freq.$}}{Percentage of times the character '$' appeared in the e-mail}
#'   \item{\code{char.freq.#}}{Percentage of times the character '#' appeared in the e-mail}
#'   \item{\code{capital.run.length.average}}{Average length of contiguous runs of capital letters in the e-mail}
#'   \item{\code{capital.run.length.longest}}{Maximum length of contiguous runs of capital letters in the e-mail}
#'   \item{\code{capital.run.length.total}}{Total number of capital letters in the e-mail}
#'}
#' @source \url{https://archive.ics.uci.edu/ml/datasets/spambase/}
"spambase"