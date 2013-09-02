

%%Contact location

tip = input('tip coordinates'); % enter lead tip coordinates relative to AC-PC [lat,AP,vert]
entry = input('entry coordinates');  % enter entry point coordinates relative to AC-PC [lat,AP,vert]
model = 3389;% input('lead model'); % 3389 or 3387
active_contact =  input('active_contact'); % negative contact
reference_contact = input('reference_contact'); % positive contact if bipolar otherwise leave blank

% length of lead in the brain
x = entry(1) - tip(1);
y = entry(2) - tip(2);
z = entry(3) - tip(3);

L = sqrt(x^ 2 + y^2 + z^2);

% distance between midpoints of adjacent contacts
if model==3387
    D = 3;
else
    D = 2;
end

% distance from tip to midpoint of active contact array;
if ~isempty(reference_contact)
    Q = 0.75 + D * (active_contact+reference_contact)/2;
else
    Q = 0.75 + D *active_contact;
end

% Coordinate of active contacts
x_contact =(Q/L)*x +tip(1);
y_contact = (Q/L)*y + tip(2); 
z_contact = (Q/L)*z +tip(3);  
contact_coordinates = [x_contact,y_contact,z_contact];

% angles
sagital_angle = atan (y_contact/z_contact); 
coronal_angle = atan (x_contact/z_contact);