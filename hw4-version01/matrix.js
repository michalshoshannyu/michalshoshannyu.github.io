const identity = () => [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];

function perspectiveMatrix(fov, aspect, near, far) {
  let f = 1.0 / Math.tan(fov / 2);
  let rangeInv = 1 / (near - far);

  return [
    f / aspect, 0, 0, 0,
    0, f, 0, 0,
    0, 0, (near + far) * rangeInv, -1,
    0, 0, near * far * rangeInv * 2, 0,
  ];
}
const translation = (x, y, z) => [
  1,
  0,
  0,
  0,
  0,
  1,
  0,
  0,
  0,
  0,
  1,
  0,
  x,
  y,
  z,
  1,
];

const rotationX = (theta) => [
  1,
  0,
  0,
  0,
  0,
  Math.cos(theta),
  Math.sin(theta),
  0,
  0,
  -Math.sin(theta),
  Math.cos(theta),
  0,
  0,
  0,
  0,
  1,
];

const rotationY = (theta) => [
  Math.cos(theta),
  0,
  -Math.sin(theta),
  0,
  0,
  1,
  0,
  0,
  Math.sin(theta),
  0,
  Math.cos(theta),
  0,
  0,
  0,
  0,
  1,
];
const rotationZ = (theta) => [
  Math.cos(theta),
  Math.sin(theta),
  0,
  0,
  -Math.sin(theta),
  Math.cos(theta),
  0,
  0,
  0,
  0,
  1,
  0,
  0,
  0,
  0,
  1,
];

const transpose = (matrix) => {
  const transposedMatrix = new Array(16);
  for (let row = 0; row < 4; row++) {
    for (let col = 0; col < 4; col++) {
      const index = col * 4 + row;
      const transposedIndex = row * 4 + col;
      transposedMatrix[transposedIndex] = matrix[index];
    }
  }
  return transposedMatrix;
};
const scale = (x, y, z) => [x, 0, 0, 0, 0, y, 0, 0, 0, 0, z, 0, 0, 0, 0, 1];

const multiply = (matrix1, matrix2) => {
  const matrix = [];
  for (let index = 0; index < 16; index++) {
    const col = Math.floor(index / 4);
    const row = index % 4;
    matrix[index] = 0;
    for (let i = 0; i < 4; i++) {
      const index1 = i * 4 + row;
      const index2 = col * 4 + i;
      matrix[index] += matrix1[index1] * matrix2[index2];
    }
  }
  return matrix;
};

