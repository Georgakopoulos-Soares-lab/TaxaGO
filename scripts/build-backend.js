const { execSync } = require('child_process');
const os = require('os');
const path = require('path');
const fs = require('fs');

const repoRoot = path.join(__dirname, '..');

const targetMap = {
  win32: 'x86_64-pc-windows-msvc',
  darwin: os.arch() === 'arm64' ? 'aarch64-apple-darwin' : 'x86_64-apple-darwin',
  linux: 'x86_64-unknown-linux-gnu'
};

const platform = process.platform;
const target = targetMap[platform];
if (!target) {
  console.error(`Unsupported platform ${platform}`);
  process.exit(1);
}

console.log('[build-backend] compiling taxago binaries for', target);

// Build all binaries
const binaries = ['taxago', 'semantic-similarity', 'common-ancestors'];
binaries.forEach(binary => {
  console.log(`[build-backend] building ${binary}...`);
  execSync(`cargo build --release --bin ${binary} --target ${target}`, { stdio: 'inherit', cwd: repoRoot });
});

// Copy all binaries to electron folder
binaries.forEach(binary => {
  const binName = platform === 'win32' ? `${binary}.exe` : binary;
  const builtPath = path.join(repoRoot, 'target', target, 'release', binName);
  const destPath = path.join(repoRoot, 'electron', binName);
  
  if (fs.existsSync(builtPath)) {
    fs.copyFileSync(builtPath, destPath);
    console.log(`[build-backend] copied ${binName} to electron folder`);
  } else {
    console.error(`[build-backend] ERROR: ${binName} not found at ${builtPath}`);
    process.exit(1);
  }
});

console.log('[build-backend] all binaries built and copied successfully');
