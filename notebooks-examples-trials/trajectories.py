from Magics.macro import *



def trajectories(title_size=0.4, legend_size=0.3):
	toplot = []
	#settings of geographical area 
	area = mmap(subpage_map_projection="cylindrical",
			subpage_lower_left_longitude=-110.,
			subpage_lower_left_latitude=20.,
			subpage_upper_right_longitude=-30.,
			subpage_upper_right_latitude=70.,
		)
	toplot.append(area)

	#settings of the coastlines attributes 
	background_coast = mcoast(
	  map_coastline_land_shade = "on",
	  map_coastline_land_shade_colour = "cream",
	  map_grid_line_style = "dash",
	  map_grid_colour = "grey",
	  map_label = "on",
	  map_coastline_colour = "grey"
	)
	toplot.append(background_coast)
	#settings of the coastlines attributes 
	foreground_coast = mcoast(
	  map_grid_line_style = "dash",
	  map_grid_colour = "grey",
	  map_label_colour = "charcoal",
	  map_coastline_colour = "grey"
	)

	  
	colours = [ "red", "blue", "green", "yellow", "purple", "orange", "cyan", "black" ]
	
	colour = 0;
	for i in range(1, 52):
			
		file = "trajectory_%02d.csv" % i
		
		data = mtable(table_filename = file,
				table_variable_identifier_type='index',
				table_latitude_variable = "1",
				table_longitude_variable = "2",
				table_value_variable = "6",
				table_header_row = 0,
				 )

		line=msymb(symbol_type='marker',
				symbol_marker_index = 28,
				symbol_colour = colours[colour],
				symbol_height = 0.20,
				symbol_text_font_colour = "black",		
				symbol_connect_line ='on'
				)
		colour += 1
		if colour == len(colours) :
			colour = 0
		toplot.append(data)
		toplot.append(line)
		


	
	title = mtext(
		text_lines= ["<font colour='navy'> Sandi possible trajectories </font>"],
		text_justification= 'centre',
		text_font_size= 1.,
		text_font_style= 'bold',
		text_box_blanking='on',
		text_border='on',
		text_mode='positional',
		text_box_x_position = 3.,
		text_box_y_position = 12.,
		text_box_x_length = 13.,
		text_box_y_length = 1.5,
		text_border_colour='charcoal')
	
	
	toplot.append(foreground_coast)
	toplot.append(title)
	return toplot
	
if __name__ == "__main__":
    #settings of the png output 
	output = output(
			output_formats = ['png'],
  			output_name = "trajectories",
    		output_name_first_page_number = "off"
    )

	plot(output, trajectories(title_size = 0.8, legend_size=0.4))
    
