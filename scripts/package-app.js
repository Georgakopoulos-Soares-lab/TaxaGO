#!/usr/bin/env node

const { execSync } = require('child_process');
const path = require('path');
const fs = require('fs');

const repoRoot = path.join(__dirname, '..');
const electronDir = path.join(repoRoot, 'electron');

console.log('🚀 TaxaGO Packaging Script');
console.log('==========================');

// Check if we're in the right directory
if (!fs.existsSync(path.join(repoRoot, 'Cargo.toml'))) {
  console.error('❌ Error: Cargo.toml not found. Please run this script from the project root.');
  process.exit(1);
}

// Check if electron directory exists
if (!fs.existsSync(electronDir)) {
  console.error('❌ Error: electron directory not found.');
  process.exit(1);
}

// Step 1: Generate custom icons
console.log('\n🎨 Step 1: Generating custom TaxaGO icons...');
try {
  execSync('node scripts/generate-icons.js', { stdio: 'inherit', cwd: repoRoot });
  console.log('✅ Custom icons generated successfully');
} catch (error) {
  console.error('❌ Failed to generate icons:', error.message);
  process.exit(1);
}

// Step 2: Install dependencies
console.log('\n📦 Step 2: Installing Electron dependencies...');
try {
  execSync('npm install', { stdio: 'inherit', cwd: electronDir });
  console.log('✅ Dependencies installed successfully');
} catch (error) {
  console.error('❌ Failed to install dependencies:', error.message);
  process.exit(1);
}

// Step 3: Build backend (this will also copy binaries)
console.log('\n🔨 Step 3: Building Rust backend...');
try {
  execSync('node ../scripts/build-backend.js', { stdio: 'inherit', cwd: electronDir });
  console.log('✅ Backend built successfully');
} catch (error) {
  console.error('❌ Failed to build backend:', error.message);
  process.exit(1);
}

// Step 4: Package the application
console.log('\n📱 Step 4: Packaging application...');
try {
  execSync('npm run dist', { stdio: 'inherit', cwd: electronDir });
  console.log('✅ Application packaged successfully');
} catch (error) {
  console.error('❌ Failed to package application:', error.message);
  process.exit(1);
}

// Step 5: Show results
console.log('\n🎉 Packaging completed successfully!');
console.log('📁 Check the electron/dist/ directory for your packaged application.');

// List the output files
const distDir = path.join(electronDir, 'dist');
if (fs.existsSync(distDir)) {
  console.log('\n📋 Generated files:');
  const files = fs.readdirSync(distDir);
  files.forEach(file => {
    const filePath = path.join(distDir, file);
    const stats = fs.statSync(filePath);
    const size = (stats.size / 1024 / 1024).toFixed(2);
    console.log(`   📄 ${file} (${size} MB)`);
  });
}

console.log('\n✨ You can now distribute your TaxaGO application!');
