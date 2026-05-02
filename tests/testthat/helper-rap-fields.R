make_rap_fields <- function() {
  data.frame(
    field_name = c(
      "participant.eid",
      "participant.p31",
      "participant.p53_i0",
      "participant.p21022",
      "participant.p21001_i0",
      "participant.p4080_i0_a0",
      "participant.p4080_i0_a1",
      "participant.p93_i0_a0",
      "participant.p20002_i0_a0"
    ),
    title = c(
      "Participant ID",
      "Sex",
      "Date of attending assessment centre | Instance 0",
      "Age at recruitment",
      "Body mass index | Instance 0",
      "Systolic blood pressure, automated reading | Instance 0 | Array 0",
      "Systolic blood pressure, automated reading | Instance 0 | Array 1",
      "Systolic blood pressure, manual reading | Instance 0 | Array 0",
      "Non-cancer illness code, self-reported | Instance 0 | Array 0"
    ),
    stringsAsFactors = FALSE
  )
}
