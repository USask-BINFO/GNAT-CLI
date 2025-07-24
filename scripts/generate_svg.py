#! /usr/bin/env python

# DEPRECATED. Was used to draw representations of GN plots based on avg orientation.

import svgwrite

def create_forward_marker(dwg, x, y, width, height, colour):
    marker = dwg.marker(
        insert=(0, 1),
        size=(2, 2),
        orient='auto'
    )
    # marker = dwg.marker(insert=(0, size), size=(size, size), orient='auto')
    marker.add(dwg.path(d="M 0 0 L 1 1 L 0 2 z", fill=colour))

    dwg.defs.add(marker)

    return marker

def create_reverse_marker(dwg, size, colour):
    marker = dwg.marker(insert=(1, 1), size=(2, 2), orient='auto')
    marker.add(dwg.path(d="M 1 0 L 0 1 L 1 2 z", fill=colour))

    dwg.defs.add(marker)

    return marker

def draw_line_with_arrow(dwg, start, length, colour, width, arrow_width, arrow_height, reverse=False):
    """
    Draws a line with an arrow at one end. Starts at start, and goes right for length.
    Width is line width. Arrow size is the size of the arrowhead. If reverse=True, the arrow
    will point left instead of right.
    """
    start_x = start[0]
    start_y = start[1]
    end_x = start_x + length - arrow_width/2
    end_y = start_y


    if reverse:
        start_x = start_x + arrow_width/2
        end_x = end_x + arrow_width/2
        line = dwg.line(start=(start_x, start_y), end=(end_x, end_y), stroke=colour, stroke_width=width)
        arrow = create_reverse_marker(dwg, arrow_width, colour)
        line['marker-start'] = arrow.get_funciri()
    else:
        line = dwg.line(start=start, end=(end_x, end_y), stroke=colour, stroke_width=width)
        arrow = create_forward_marker(dwg, end_x, end_y, arrow_width, arrow_height, colour)
        line['marker-end'] = arrow.get_funciri()
    
    dwg.add(line)

def convert_orientations_ones(orientation):
    # Convert orientation from [1, 1, -1, 1, 0, 1, 1, -1, -1, -1, 1] to [[2, -1],[4],[-3, 1]]
    orientation = [int(x) for x in orientation if x != 'None'] # Trim the Nones off
    orientation_summary = []
    prev = orientation[0]
    current = None
    current_block = [prev]
    i = 1
    while i < len(orientation):
        current = orientation[i]
        if current == prev:
            current_block.append(current)
            prev = current
        else:
            if current == 0:
                if prev == 1:
                    current_block.append(current)
                else:
                    orientation_summary.append(current_block)
                    current_block = [current]
                    prev = 1
            else:
                orientation_summary.append(current_block)
                current_block = [current]
                prev = current
        i += 1

        if i == len(orientation):
            orientation_summary.append(current_block)

    before, match, after = [], [], []
    for i in range(len(orientation_summary)):
        if 0 in orientation_summary[i]:
            match += [len(orientation_summary[i])]
        else:
            if len(match) == 0:
                before += [len(orientation_summary[i]) * -1 if orientation_summary[i][0] < 0 else len(orientation_summary[i])]
            else:
                after += [len(orientation_summary[i]) * -1 if orientation_summary[i][0] < 0 else len(orientation_summary[i])]

    orientation_summary = [before, match, after]

    return orientation_summary

def create_svg_from_orientation(orientation, filename):
    # Convert orientation from [1, 1, -1, 1, 0, 1, 1, -1, -1, -1, 1] to [[2, -1],[4],[-3, 1]]
    before, match, after = convert_orientations_ones(orientation)

    total_length = sum(10*abs(val) for val in before) + sum(10*abs(val) for val in match) + sum(10*abs(val) for val in after)

    # Define SVG drawing parameters
    image_width = total_length + 10
    image_height = 50
    stroke = 5
    arrow_width = 10
    arrow_height = stroke * 1.5
    forward_arrow_colour = 'green'
    reverse_arrow_colour = 'red'
    match_arrow_colour = 'gold'
    y = image_height / 2

    # Begin the drawing
    dwg = svgwrite.Drawing(filename, size=(image_width, image_height))

    # Draw before
    x = 10
    for i in range(len(before)):
        abs_before = abs(before[i])
        if before[i] < 0:
            draw_line_with_arrow(dwg, (x, y), abs_before*10, reverse_arrow_colour, stroke, arrow_width, arrow_height, reverse=True)
        else:
            draw_line_with_arrow(dwg, (x, y), abs_before*10, forward_arrow_colour, stroke, arrow_width, arrow_height)
        x += abs_before *10

    # Draw match
    draw_line_with_arrow(dwg, (x, y), match[0]*10, match_arrow_colour, stroke, arrow_width, arrow_height)
    x += abs(match[0])*10

    # Draw after
    for i in range(len(after)):
        abs_after = abs(after[i])
        if after[i] < 0:
            draw_line_with_arrow(dwg, (x, y), abs_after*10, reverse_arrow_colour, stroke, arrow_width, arrow_height, reverse=True)
        else:
            draw_line_with_arrow(dwg, (x, y), abs_after*10, forward_arrow_colour, stroke, arrow_width, arrow_height)
        x += abs_after*10

    dwg.save()

