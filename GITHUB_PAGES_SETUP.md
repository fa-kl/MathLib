# Setting Up GitHub Pages for MathLib Documentation

This guide will walk you through setting up automatic documentation deployment to GitHub Pages.

## Step 1: Enable GitHub Pages

1. Go to your repository: https://github.com/fa-kl/MathLib
2. Click on **Settings** (top navigation bar)
3. Scroll down to **Pages** in the left sidebar
4. Under **Source**, select:
   - Source: **GitHub Actions**
   
That's it! The GitHub Actions workflow is already configured in `.github/workflows/deploy-docs.yml`

## Step 2: Push Your Changes

```bash
# Add all the new files
git add .

# Commit with a descriptive message
git commit -m "Release v1.0.0: Add versioning, documentation, and GitHub Pages setup"

# Push to main branch
git push origin main
```

## Step 3: Trigger Documentation Build

The documentation will build automatically when you push to `main`. You can also:

1. Go to **Actions** tab in your GitHub repository
2. Click on **Deploy Documentation** workflow
3. Click **Run workflow** button
4. Select the `main` branch
5. Click **Run workflow**

## Step 4: Access Your Documentation

After the workflow completes (usually 2-3 minutes):

**Your documentation will be available at:**
### 🌐 https://fa-kl.github.io/MathLib

## How It Works

1. **Doxyfile**: Configures Doxygen to parse your code comments
2. **GitHub Actions Workflow**: Automatically runs Doxygen and deploys to GitHub Pages
3. **On Every Push**: Documentation updates automatically when you push to main

## Updating Documentation

Simply push changes to your code comments or README, and the documentation will rebuild automatically!

## Local Documentation Generation

To generate documentation locally:

```bash
# Install Doxygen (if not already installed)
sudo apt-get install doxygen graphviz

# Generate documentation
doxygen Doxyfile

# Open in browser
xdg-open docs/html/index.html
```

## Troubleshooting

### Documentation not showing up?
1. Check the Actions tab for any errors
2. Make sure GitHub Pages is enabled in Settings
3. Wait a few minutes for DNS propagation

### Build failing?
1. Check the workflow logs in the Actions tab
2. Ensure Doxyfile is valid
3. Verify all code comments use proper Doxygen syntax

---

**Note**: The first deployment may take 5-10 minutes. Subsequent updates are faster.
