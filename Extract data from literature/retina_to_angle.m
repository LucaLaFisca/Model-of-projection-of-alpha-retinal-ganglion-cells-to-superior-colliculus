function y = retina_to_angle(size_on_retina)
eye_length = 3.37;
y = atan(size_on_retina/eye_length)*180/pi;