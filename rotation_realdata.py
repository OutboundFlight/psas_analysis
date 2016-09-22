"""
	PSAS Orientation Algorithm

	Originally developed by Aaron Baker
	Translated to python by William Harrington
"""
# panda, panda
import pandas as pd

# for matrix and math stuff
import numpy as np

# quaternion stuff
from thLib import quat

# read in csv file
data = pd.read_csv('ADIS.csv')

# extract timestamps
time = data[' [1]Timestamp']

# umblicial disconnect event
t_0 = 117853569585227
time = np.subtract(time, t_0)

# convert from ns to s
time = np.divide(time, 1e9)

# make numpy array (plays nicer)
time = np.array(time[3650:28125])

# create vector from data, eliminating samples from before takeoff and after apogee
# Note: key has space because of our csv formatting mess
gyro_data = zip(np.deg2rad(data[' [3]Gyro_X'][3650:28125]),
                np.deg2rad(data[' [4]Gyro_Y'][3650:28125]),
                np.deg2rad(data[' [5]Gyro_Z'][3650:28125]))
gyro_data = np.transpose(gyro_data)

steps = len(time)
tstep = np.mean(np.diff(time))

orientations = np.zeros((4, steps))
angles = np.zeros((3, steps))
dcm = np.zeros((3,3,steps))

orientations[:,0] = np.transpose([1,0,0,0])
angles[:,0] = np.deg2rad(quat.quat2deg(orientations[:,0]))

""" TODO: Finish translation """
"""
dcm(1:3,1:3,1) = quat2dcm(orientations(1:4,1)');


%from the book "Strapdown Intertial Navigation Technology", section 13.7
for i = 2:steps

    %Multiply the rates from the gyro by the time step length to get a
    %rotation amount ("Zero-order hold assumption")
    r = tstep*rates(1:3,i);

    %Sigma is the magnitude of the vector of roll rates from the rate gyro,
    %times the length of the time step
    halfsigmasq = 0.25*(r(1)*r(1) + r(2)*r(2) + r(3)*r(3));
    term2 = halfsigmasq*halfsigmasq;
    
    %modified taylor series expansions of sine and cosine
    sac = 1 - halfsigmasq/2 + term2/24;
    sas = 0.5*(1 - halfsigmasq/6 + term2/120);    

    %the update quaternion as derived in the text; derivation is very
    %mysterious
    q = [sac, r(1)*sas, r(2)*sas, r(3)*sas];

    %the updated orientatoin is the previous orientation times the update
    %quaternion
    o = quatmultiply(orientations(1:4, i-1)', q);

    %renormalize the orientation at each time step by dividing by the
    %length
    o = o/norm(o)';
    orientations(1:4, i) = o(1:4)';
    
    %calculate the euler angle and DCM equivalents
    [angles(1,i), angles(2,i), angles(3,i)]  = quat2angle(orientations(1:4,i)');
    angles(1:3,i) = rad2deg*angles(1:3,i);
    dcm(1:3,1:3,i) = quat2dcm(orientations(1:4,i)');
    
end

%Grab the magnetic field vector directly from the data by normalizing the
%xyz axis magnetometer measurements
magn = zeros(3, steps);
magrot = zeros(3,steps);
for i = 1:steps
    m = [data(i, 10), data(i,11), data(i,12)]';
    magn(:, i) = m/norm(m);
    magrot(:,i) = quatrotate(quatinv(orientations(1:4, i)'), magn(:,i)');
end

%reshape matrices for quiver plots 
x1 = reshape(dcm(1,1,:), 1, steps);
y1 = reshape(dcm(2,1,:), 1, steps);
z1 = reshape(dcm(3,1,:), 1, steps);

x2 = reshape(dcm(1,2,:), 1, steps);
y2 = reshape(dcm(2,2,:), 1, steps);
z2 = reshape(dcm(3,2,:), 1, steps);

x3 = reshape(dcm(1,3,:), 1, steps);
y3 = reshape(dcm(2,3,:), 1, steps);
z3 = reshape(dcm(3,3,:), 1, steps);

%number of data points to skip between frames of the movie
skip = 100;

F(length(1:skip:steps)) = struct('cdata', [], 'colormap', []);
index = 1;

%Plot the DCM and magnetometer vectors as individual arrows on a quiver
%plot and save frames to make a movie
for i = 1:skip:steps
    quiver3(0,0,0, x1(i), y1(i), z1(i), 0);
    hold on
    title(sprintf('Time: %g sec\t Accel: %g sec', data(i, 2), norm([data(i,7), data(i, 8), data(i,9)])))
    quiver3(0,0,0, x2(i), y2(i), z2(i), 0);
    quiver3(0,0,0, x3(i), y3(i), z3(i), 0);
    quiver3(0,0,0, magn(1,i), magn(2,i), magn(3,i), 0);
    %quiver3(0,0,0, magrot(1,i), magrot(2,i), magrot(3,i), 0);
    hold off
    xlim([-1,1]);
    ylim([-1,1]);
    zlim([-1,1]);
    legend('Local X', 'Local Y', 'Local Z', 'Magnetic Field direction')
    
    F(index) = getframe;
    index = index + 1;
end

%plot3(x1,y1,z1);
%grid on
%hold on
%plot3(x2,y2,z2);
%plot3(x3,y3,z3);
%hold off

%angdist=zeros(1,steps);
%dotprod = zeros(1,steps);
%for i = 1:steps
%    v = [x2(i) y2(i) z2(i)];
%    dotprod(i) = dot(magn(:,i)', v);
%    angdist(i) = acosd(dotprod(i));
%end
"""
