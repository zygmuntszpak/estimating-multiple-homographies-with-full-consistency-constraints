e1 = [1,0,0]';
e2 = [0,1,0]';
e3 = [0,0,1]';

E1 = [e2*e1' + e3*e2', e3*e3' - e1*e1', -e1*e2' - e2*e3']
E2 = [e1*e3' - e2*e2' + e3*e1']
rand( 'seed', 0 );
randn( 'seed', 10 );

A = rand(3,3);
B = rand(3,3);
C = -E1*(kron(eye(3),E2))*kron(A,B)

x = rand(3,1);
y = rand(3,1);

first = cross(A*x,B*y);
second = C*(kron(x,y))

first - second

D = [1 0 0 0 0 0;
     0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 1 0 0;
     0 1 0 0 0 0;
     0 0 0 0 0 1;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     0 0 1 0 0 0];
 
 l = [1 2 3]';
 l = rand(3,1);
 
 one = kron(l,l)
 two = D*[l(1)^2 l(2)^2 l(3)^2 l(1)*l(2) l(1)*l(3) l(2)*l(3)]'
 
 one - two
 
 test1 = l*l'
 test2 = create_matrix_N([l(1)^2 l(2)^2 l(3)^2 l(1)*l(2) l(1)*l(3) l(2)*l(3)]')
 
 test1 - test2
 
 