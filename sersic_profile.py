def circle_intensity(array, centre , radius):
    max_y, min_y = centre[1] + radius, centre[1] - radius
    max_x, min_x = centre[0] + radius, centre[0] - radius
    length = len(array)
    cumulative_sum = 0
    for i in range(min_y,max_y + 1):
        for j in range(min_x, max_x + 1):
            if ( (i - centre[1]) **2  + (j - centre[0])**2 ) <= (radius * radius):
                cumulative_sum += array[i][j]                
    return cumulative_sum 


def half_light_radius(array, centre, full_radius, total_intensity):
    half_intensity = total_intensity * 0.5
    previous_intensity = 0
    for i in range(1,full_radius):
        current_intensity = circle_intensity(array, centre, i)
        if current_intensity >= half_intensity:
            correction = (current_intensity - half_intensity) / (current_intensity - previous_intensity)    # Linear correction
            return (i, i - correction)          # returns integer radius and linearly corrected radius 
        previous_intensity = current_intensity
