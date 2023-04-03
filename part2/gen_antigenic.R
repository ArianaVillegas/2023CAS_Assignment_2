library("Racmacs")
path_to_titer_file <- "part2/titer_escape_b10_d2.csv"    # Get the file path
titer_table        <- read.titerTable(path_to_titer_file) # Read into table
map <- make.acmap(                                        # Generate map
	titer_table             = titer_table,
	number_of_dimensions    = 2,
	number_of_optimizations = 500,
	minimum_column_basis    = "none"
)

pdf("antigenic_map_escape_b10_d2.pdf") # Open saving device
plot(map)                               # Create the plot
dev.off()                               # Close the device
