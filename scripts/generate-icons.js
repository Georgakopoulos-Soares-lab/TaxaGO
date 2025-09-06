#!/usr/bin/env node

const { execSync } = require('child_process');
const path = require('path');
const fs = require('fs');

const repoRoot = path.join(__dirname, '..');
const assetsDir = path.join(repoRoot, 'electron', 'assets');
const logoPath = path.join(repoRoot, 'gorilla.svg');

console.log('🎨 TaxaGO Icon Generation Script');
console.log('================================');

// Check if logo exists
if (!fs.existsSync(logoPath)) {
  console.error('❌ Error: gorilla.svg not found in project root');
  process.exit(1);
}

// Check if assets directory exists
if (!fs.existsSync(assetsDir)) {
  console.error('❌ Error: electron/assets directory not found');
  process.exit(1);
}

console.log('📁 Assets directory:', assetsDir);
console.log('🎯 Gorilla logo file:', logoPath);

// Function to check if a command exists
function commandExists(command) {
  try {
    execSync(`which ${command}`, { stdio: 'ignore' });
    return true;
  } catch {
    return false;
  }
}

// Check for required tools
const requiredTools = ['rsvg-convert', 'convert'];
const missingTools = [];

requiredTools.forEach(tool => {
  if (!commandExists(tool)) {
    missingTools.push(tool);
  }
});

if (missingTools.length > 0) {
  console.error('❌ Missing required tools:', missingTools.join(', '));
  console.log('\n📦 Please install the missing tools:');
  console.log('   macOS: brew install librsvg imagemagick');
  console.log('   Ubuntu/Debian: sudo apt-get install librsvg2-bin imagemagick');
  console.log('   Windows: Install ImageMagick and librsvg');
  process.exit(1);
}

console.log('✅ All required tools found');

// Generate PNG icons from SVG with rounded corners
console.log('\n🖼️  Generating PNG icons with rounded corners...');

const sizes = [16, 32, 64, 128, 256, 512];
const tempIconsDir = path.join(assetsDir, 'temp_icons');

// Create temp directory if it doesn't exist
if (!fs.existsSync(tempIconsDir)) {
  fs.mkdirSync(tempIconsDir, { recursive: true });
}

sizes.forEach(size => {
  const outputPath = path.join(tempIconsDir, `icon_${size}.png`);
  const tempPath = path.join(tempIconsDir, `temp_${size}.png`);
  console.log(`   Generating ${size}x${size} icon with rounded corners...`);
  
      try {
      // Make the overall icon smaller by using 85% of the target size for the background
      const backgroundSize = Math.floor(size * 0.85);
      const backgroundOffset = Math.floor((size - backgroundSize) / 2);
      
      // Generate the gorilla at 90% of the background size (so it's 76.5% of the total size)
      const iconSize = Math.floor(backgroundSize * 0.9);
      const iconOffset = backgroundOffset + Math.floor((backgroundSize - iconSize) / 2);
      
      // Create a temporary file with the icon
      execSync(`rsvg-convert -w ${iconSize} -h ${iconSize} "${logoPath}" -o "${tempPath}"`, { stdio: 'inherit' });
      
      // Create a smaller rounded rectangle background with the same color as the SVG and composite the icon on top
      const radius = Math.floor(backgroundSize * 0.15); // 15% of background size for corner radius
      
      execSync(`convert -size ${size}x${size} xc:transparent -fill "#f8fafc" -draw "roundrectangle ${backgroundOffset},${backgroundOffset} ${backgroundOffset + backgroundSize},${backgroundOffset + backgroundSize} ${radius},${radius}" "${tempPath}" -geometry +${iconOffset}+${iconOffset} -composite "${outputPath}"`, { stdio: 'inherit' });
    
    // Clean up temp file
    fs.unlinkSync(tempPath);
    
    console.log(`   ✅ Generated ${outputPath}`);
  } catch (error) {
    console.error(`   ❌ Failed to generate ${size}x${size} icon:`, error.message);
  }
});

// Create main icon.png (512x512) with rounded corners
const mainIconPath = path.join(assetsDir, 'icon.png');
console.log(`\n📄 Creating main icon.png with rounded corners...`);
try {
  // Make the overall icon smaller by using 85% of the target size for the background
  const backgroundSize = Math.floor(512 * 0.85);
  const backgroundOffset = Math.floor((512 - backgroundSize) / 2);
  
  // Generate the gorilla at 90% of the background size (so it's 76.5% of the total size)
  const iconSize = Math.floor(backgroundSize * 0.9);
  const iconOffset = backgroundOffset + Math.floor((backgroundSize - iconSize) / 2);
  
  const radius = Math.floor(backgroundSize * 0.15);
  const tempPath = path.join(assetsDir, 'temp_main.png');
  
  // Create a temporary file with the icon
  execSync(`rsvg-convert -w ${iconSize} -h ${iconSize} "${logoPath}" -o "${tempPath}"`, { stdio: 'inherit' });
  
  // Create a smaller rounded rectangle background with the same color as the SVG and composite the icon on top
  execSync(`convert -size 512x512 xc:transparent -fill "#f8fafc" -draw "roundrectangle ${backgroundOffset},${backgroundOffset} ${backgroundOffset + backgroundSize},${backgroundOffset + backgroundSize} ${radius},${radius}" "${tempPath}" -geometry +${iconOffset}+${iconOffset} -composite "${mainIconPath}"`, { stdio: 'inherit' });
  
  // Clean up temp file
  fs.unlinkSync(tempPath);
  
  console.log('✅ Created main icon.png with rounded corners');
} catch (error) {
  console.error('❌ Failed to create main icon.png:', error.message);
}

// Create Windows ICO file
console.log('\n🪟 Creating Windows ICO file...');
try {
  const iconFiles = sizes.map(size => path.join(tempIconsDir, `icon_${size}.png`)).join(' ');
  execSync(`convert ${iconFiles} "${path.join(assetsDir, 'icon.ico')}"`, { stdio: 'inherit' });
  console.log('✅ Created icon.ico');
} catch (error) {
  console.error('❌ Failed to create ICO file:', error.message);
}

// Create macOS ICNS file
console.log('\n🍎 Creating macOS ICNS file...');
try {
  // Create iconset directory
  const iconsetDir = path.join(assetsDir, 'icon.iconset');
  if (!fs.existsSync(iconsetDir)) {
    fs.mkdirSync(iconsetDir, { recursive: true });
  }

  // Copy icons to iconset with proper naming
  const iconsetMap = {
    16: 'icon_16x16.png',
    32: 'icon_16x16@2x.png',
    32: 'icon_32x32.png',
    64: 'icon_32x32@2x.png',
    128: 'icon_128x128.png',
    256: 'icon_128x128@2x.png',
    256: 'icon_256x256.png',
    512: 'icon_256x256@2x.png',
    512: 'icon_512x512.png',
    1024: 'icon_512x512@2x.png'
  };

  Object.entries(iconsetMap).forEach(([size, filename]) => {
    const sourcePath = path.join(tempIconsDir, `icon_${size}.png`);
    const destPath = path.join(iconsetDir, filename);
    if (fs.existsSync(sourcePath)) {
      fs.copyFileSync(sourcePath, destPath);
    }
  });

  // Create ICNS file
  execSync(`iconutil -c icns "${iconsetDir}" -o "${path.join(assetsDir, 'icon.icns')}"`, { stdio: 'inherit' });
  
  // Clean up iconset directory
  fs.rmSync(iconsetDir, { recursive: true, force: true });
  
  console.log('✅ Created icon.icns');
} catch (error) {
  console.error('❌ Failed to create ICNS file:', error.message);
}

console.log('\n🎉 Icon generation completed!');
console.log('📁 Generated files:');
console.log('   - icon.png (512x512)');
console.log('   - icon.ico (Windows)');
console.log('   - icon.icns (macOS)');
console.log('   - temp_icons/ (various sizes)');

console.log('\n✨ Your TaxaGO app will now use custom icons instead of the default Electron icon!');
