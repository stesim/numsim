function visualize(filename)
	file = fopen(filename);
	gridSize = fscanf(file, '%d %d');
	gridSize = gridSize(1:2);
	entries = fscanf(file, '%f %f %f %f %f', [5,Inf]);
	x = reshape(entries(1,:),[gridSize(2),gridSize(1)]).';
	y = reshape(entries(2,:),[gridSize(2),gridSize(1)]).';
	u = reshape(entries(3,:),[gridSize(2),gridSize(1)]).';
	v = reshape(entries(4,:),[gridSize(2),gridSize(1)]).';
	p = reshape(entries(5,:),[gridSize(2),gridSize(1)]).';

	subplot(2,2,1);
	surf(x,y,u);
	title('u');
	xlabel('x');
	xlabel('y');

	subplot(2,2,2);
	surf(x,y,v);
	title('v');
	xlabel('x');
	xlabel('y');

	subplot(2,2,3);
	surf(x,y,p);
	title('p');
	xlabel('x');
	xlabel('y');
end
