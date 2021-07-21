library(track2KBA)

## Test whether csv downloaded from Movebank converts correctly
expect_true( "ID" %in% colnames(move2KBA(
  filename="movebank_testdata.csv")$data),
  info = "ID column conversion from movebank file"
)

expect_true( "DateTime" %in% colnames(move2KBA(
  filename="movebank_testdata.csv")$data),
  info = "DateTime column conversion from movebank file"
)

expect_true( "Latitude" %in% colnames(move2KBA(
  filename="movebank_testdata.csv")$data),
  info = "DateTime column conversion from movebank file"
)

expect_true( "Longitude" %in% colnames(move2KBA(
  filename="movebank_testdata.csv")$data),
  info = "DateTime column conversion from movebank file"
)

expect_message(move2KBA(
  filename="movebank_testdata.csv")$data,
  info = "Colony location message"
)
