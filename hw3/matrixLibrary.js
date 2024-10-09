function Matrix() {
    this.identity = () =>
      (value = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]);
    
    this.translate = (x, y, z) =>
      (value = multiply(value, [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, x, y, z, 1]));
  
    this.rotateX = (theta) => {
      let c = Math.cos(theta);
      let s = Math.sin(theta);
      let xRotation = [1, 0, 0, 0, 0, c, s, 0, 0, -s, c, 0, 0, 0, 0, 1];
      value = multiply(value, xRotation);
    };
  
    this.rotateY = (theta) => {
      let c = Math.cos(theta);
      let s = Math.sin(theta);
      let yRotation = [c, 0, -s, 0, 0, 1, 0, 0, s, 0, c, 0, 0, 0, 0, 1];
      value = multiply(value, yRotation);
    };
  
    this.rotateZ = (theta) => {
      let c = Math.cos(theta);
      let s = Math.sin(theta);
      let zRotation = [c, s, 0, 0, -s, c, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];
      value = multiply(value, zRotation);
    };
    this.scale = (x, y, z) => {
      value = multiply(value, [x, 0, 0, 0, 0, y, 0, 0, 0, 0, z, 0, 0, 0, 0, 1]);
    };
    this.perspective = (x, y, z) => {
      value = multiply(value, [1, 0, 0, x, 0, 1, 0, y, 0, 0, 1, z, 0, 0, 0, 1]);
    };
  
    this.get = () => value;
    this.set = (v) => (value = v);
  
    this.transform = (vector) => {
      let vecTransform = [];
      for (let i = 0; i < 4; i++) {
        vecTransform[i] =
          value[i * 4 + 0] * vector[0] +
          value[i * 4 + 1] * vector[1] +
          value[i * 4 + 2] * vector[2] +
          value[i * 4 + 3] * vector[3];
      }
      return vecTransform;
    };
  
    let value = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];
  
    let multiply = (matrix1, matrix2) => {
      let matrix3 = [];
      for (let i = 0; i < 4; i++) {
        for (let j = 0; j < 4; j++) {
          matrix3[i * 4 + j] =
            matrix1[i * 4] * matrix2[j] +
            matrix1[i * 4 + 1] * matrix2[j + 4] +
            matrix1[i * 4 + 2] * matrix2[j + 8] +
            matrix1[i * 4 + 3] * matrix2[j + 12];
        }
      }
      return matrix3;
    };
  }
  